/*
 * singleProbeChaos.cpp
 *
 *  Created on: Feb 14, 2023
 *      Author: norms
 */

#include "../headers/PlanetAndProbe.h"
#include "../headers/vector2.h"
#include "../headers/circleGenerator.h"
#include "../headers/singleProbeChaos.h"
#include "../headers/toFile.h"

// Unit testing
#include "../headers/unitTests.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>

// Main function, extract two lyapunov exponents given a single probe with perturbations
double *singleProbeChaos(Planet *earth, Planet sun, Probe *jw, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, double RADIUS)
{
	// Write the circle of points to file
	writeCircleToFile(NUMPERTURBED, RADIUS);

	// Initialize perturbed probes about james webb probe
	initializePerturbations(jw, NUMPERTURBED);

	// Allocate h2 probes
	jw->h2 = (Probe *)malloc(2 * sizeof(Probe));

	// Copy info from James Webb into h2 probes
	for (int i = 0; i < 2; ++i)
	{
		copyProbe(&jw->h2[i], *jw);
	}

	// Initialize h2 orthonormal vectors:
	jw->h2[0].pos.x += RADIUS;
	jw->h2[1].pos.y += RADIUS;

	// Storage for h1 (max) and h2 (min) lyapunov exponents
	double *lyapunovs = (double *)malloc(2 * sizeof(double));

	// Sum of log of relative separation every timestep for each perturbation: ORBITAL SEPARATION METHOD FOR CALCULATING LYAPUNOVS
	// Max lyapunov will be extracted from this at the end

	// Allocate the space
	double *LfSums = (double *)malloc(NUMPERTURBED * sizeof(double));

	// Initialize all sums to zero
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		LfSums[i] = 0;
	}

	// h2 lyapunov will be extracted from this sum at the end
	double areaSum = 0;
	double initialArea = RADIUS * RADIUS;

	// Store data: Store all data of probe and perturbations every ITERATIONS step: Useful for plotting.
	// Number of data buckets needed
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;

	// Position Data for earth and james webb
	vector2 *earthData = (vector2 *)malloc(n * sizeof(vector2));
	vector2 *jwData = (vector2 *)malloc(n * sizeof(vector2));

	// Data for perturbations: 2d array
	vector2 **perturbedData = (vector2 **)malloc(n * sizeof(vector2 *));

	// Data for h2: 2d array
	vector2 **h2Data = (vector2 **)malloc(n * sizeof(vector2 *));

	// Allocate memory for all 2d data arrays
	for (int i = 0; i < n; ++i)
	{
		perturbedData[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		h2Data[i] = (vector2 *)malloc(2 * sizeof(vector2));
	}

	// Loop
	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{
		// Store Data in different arrays
		if (c % ITERATIONS == 0)
		{
			std::cout << "c is: " << c << std::endl;
			// Earth
			earthData[dc] = {earth->pos.x, earth->pos.y};

			// James Webb
			jwData[dc] = {jw->pos.x, jw->pos.y};

			// Perturbed probes
			for (int i = 0; i < NUMPERTURBED; ++i)
			{
				perturbedData[dc][i] = {jw->perturbed[i].pos.x, jw->perturbed[i].pos.y};
			}

			// h2
			for (int i = 0; i < 2; ++i)
			{
				h2Data[dc][i] = {jw->h2[i].pos.x, jw->h2[i].pos.y};
			}

			// Increment the data count
			dc++;
		}

		// Numerical Integrator: RK4
		singleProbeRK4(earth, &sun, jw, NUMPERTURBED, dt);

		// Function for LfSums, Gram Schmidt2, and Renormalization of everything: NOT TESTED YET
		// lyapunovChaosStuff(jw, NUMPERTURBED, LfSums, RADIUS, &areaSum, initialArea, RADIUS);
		// lyapunovChaosStuffTest(jw, NUMPERTURBED, LfSums, RADIUS, &areaSum, initialArea, RADIUS);

		// Increment timestep count
		c++;
	}

	// Store latest data:

	// Earth
	earthData[n - 1].x = earth->pos.x; earthData[n - 1].y = earth->pos.y;

	// James webb
	jwData[n - 1].x = jw->pos.x; jwData[n - 1].y = jw->pos.y;

	// Perturbed probes
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		perturbedData[n - 1][i] = {jw->perturbed[i].pos.x, jw->perturbed[i].pos.y};
	}

	// h2
	for (int i = 0; i < 2; ++i)
	{
		h2Data[n - 1][i] = {jw->h2[i].pos.x, jw->h2[i].pos.y};
	}

	// Calculate and store lyapunov exponents: NOT TESTED
	extractLyapunovs(LfSums, areaSum, NUMPERTURBED, lyapunovs, dt * (double)totalSteps);

	// Write data to file
	h2Evolution(jwData, h2Data, n);
	//perturbedEvolution(jwData, perturbedData, n, NUMPERTURBED);
//	earthJWdataToFile(earthData, jwData, perturbedData, h2Data, n, NUMPERTURBED);

	// Free stuff

	// Free all data arrays
	for (int i = 0; i < n; ++i)
	{
		// Perturbed data
		free(perturbedData[i]);

		// h2 data
		free(h2Data[i]);
	}

	free(perturbedData);
	free(h2Data);

	// Earth and James Webb data
	free(earthData);
	free(jwData);

	// Free Lfsums
	free(LfSums);

	// Free perturbed probes of James Webb
	free(jw->perturbed);

	// Free h2 probes of James Webb
	free(jw->h2);

	// Return the calculated lyapunov exponents
	return lyapunovs;
}

// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovs(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t)
{
	// Calculate and store h1: max lyapunov

	// Find maximum sum of all perturbed probes
	double max = LfSums[0];
	for (int i = 1; i < NUMPERTURBED; ++i)
	{
		if (LfSums[i] > max)
		{
			max = LfSums[i];
		}
	}

	// The maximum lyapunov exponent is the maximum divided by the total elapsed time
	lyapunovs[0] = max / t;

	// Calculate and store the second lyapunov exponent
	lyapunovs[1] = (areaSum / t) - lyapunovs[0];
}

// Handles the timestep processes for calculating h1 and h2 lyapunov exponents
void lyapunovChaosStuff(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea, double RADIUS)
{
	// h1: Max lyapunov
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Calculate Lf of each perturbation, and renormalize the position vector
		vector2 v = {jw->perturbed[i].pos.x - jw->pos.x, jw->perturbed[i].pos.y - jw->pos.y};

		// Scale by initial distance
		double dist = mag(v) / RADIUS;
		double logDist = log(dist);

		// Add to the running sum
		LfSums[i] += logDist;

		// Renormalize the perturbed probe preserving direction
		vector2 uv = unitVector(v);

		// Scale the vector by initial distance
		uv.x *= RADIUS;
		uv.y *= RADIUS;

		jw->perturbed[i].pos.x = jw->pos.x + uv.x; jw->perturbed[i].pos.y = jw->pos.y + uv.y;

		// RADIUS > 1.0
//		 LfSums[i] += log(mag(v) / RADIUS);
//		 jw.perturbed[i].pos = {(uv.x + jw.pos.x) * RADIUS, (uv.y + jw.pos.y) * RADIUS, (uv.z + jw.pos.z) * RADIUS};
	}

	// h2 Lyapunov
	vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y};
	vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};

	// Calculate area of parallelogram: A = ||a|| ||b|| sin(ab)
	double area = mag(v1) * mag(v2) * sin(angle(v1, v2));

	double scaledArea = area / initialArea;
	// Take the log: Scale by initial area
	double logArea = log(scaledArea);

	// Add to the running sum
	(*areaSum) += logArea;
//	std::cout << std::setprecision(10) << area << std::endl;
//	std::cout << "magv1: " << mag(v1) << " magv2: " << mag(v2) << std::endl;

	// Reorthonormalize h2 vectors
	gramSchmidtRenormalizationh2(jw, NUMPERTURBED, &v1, &v2, RADIUS);
}

// Gram schmidt to re-orthonormalize unit vectors
void gramSchmidtRenormalizationh2(Probe *jw, int NUMPERTURBED, vector2 *v1, vector2 *v2, double RADIUS)
{
	// Calculate orthonormal vectors e1 and e2
	vector2 e1 = unitVector(*v1);											// First orthogonal normalized basis vector
	double v2Dote1 = dotProd(*v2, e1);										// Dot prodoct of second vector with orthobasis vector
	vector2 y2 = {v2->x - (v2Dote1) * e1.x, v2->y - (v2Dote1) * e1.y};		// Get the second orthogonal basis vector
	vector2 e2 = unitVector(y2);											// Normalize to get second orthonormal vector

	// Scale by initial distance
	e1.x *= RADIUS;
	e1.y *= RADIUS;
	e2.x *= RADIUS;
	e2.y *= RADIUS;

//	std::cout << "e1: " << e1.x << " " << e1.y << std::endl;
//	std::cout << "e2: " << e2.x << " " << e2.y << std::endl << std::endl;

	// Assign new vectors
	jw->h2[0].pos.x = jw->pos.x + e1.x; jw->h2[0].pos.y = jw->pos.y + e1.y;
	jw->h2[1].pos.x = jw->pos.x + e2.x; jw->h2[1].pos.y = jw->pos.y + e2.y;
}
// RK4 for sun earth system with single james webb with perturbations
void singleProbeRK4(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, double dt)
{
	// Initialize earth Copies: k2 - k4
	Planet k2Earth, k3Earth, k4Earth;

	// Initialize jw Copies: k2 - k4
	Probe k2jw, k3jw, k4jw;

	// Initialize kAccels and kVels: Index 0 is for earth, Index 1 is for probe, last 2 are for h2 probes

	// Earth:
	vector2 k1EarthAccel, k2EarthAccel, k3EarthAccel, k4EarthAccel, weightedEarthAccel;
	vector2 k1EarthVel, k2EarthVel, k3EarthVel, k4EarthVel, weightedEarthVel;

	// James Webb:
	vector2 k1jwAccel, k2jwAccel, k3jwAccel, k4jwAccel, weightedjwAccel;
	vector2 k1jwVel, k2jwVel, k3jwVel, k4jwVel, weightedjwVel;

	// Perturbations:
	vector2 *k1PerturbedAccels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k2PerturbedAccels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k3PerturbedAccels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k4PerturbedAccels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *weightedPerturbedAccels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));

	vector2 *k1PerturbedVels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k2PerturbedVels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k3PerturbedVels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *k4PerturbedVels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
	vector2 *weightedPerturbedVels = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));

	// h2 Probes:
	vector2 *k1h2Accels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k2h2Accels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k3h2Accels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k4h2Accels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *weightedh2Accels = (vector2 *)malloc(2 * sizeof(vector2));

	vector2 *k1h2Vels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k2h2Vels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k3h2Vels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *k4h2Vels = (vector2 *)malloc(2 * sizeof(vector2));
	vector2 *weightedh2Vels = (vector2 *)malloc(2 * sizeof(vector2));

	// Copy data into new copies

	// Earth:
	copyPlanet(&k2Earth, *earth);
	copyPlanet(&k3Earth, *earth);
	copyPlanet(&k4Earth, *earth);

	// James Webb:
	copyProbe(&k2jw, *jw);
	copyProbe(&k3jw, *jw);
	copyProbe(&k4jw, *jw);

	// Perturbed probes:

	// Allocate memory to perturbation arrays
	k2jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
	k3jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
	k4jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));

	// Copy over data to perturbed probes and h2 probes
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		copyProbe(&k2jw.perturbed[i], jw->perturbed[i]);
		copyProbe(&k3jw.perturbed[i], jw->perturbed[i]);
		copyProbe(&k4jw.perturbed[i], jw->perturbed[i]);
	}

	// h2 Probes:

	// Allocate memory for h2 probes
	k2jw.h2 = (Probe *)malloc(2 * sizeof(Probe));
	k3jw.h2 = (Probe *)malloc(2 * sizeof(Probe));
	k4jw.h2 = (Probe *)malloc(2 * sizeof(Probe));

	// Copy over data to h2 probes
	for (int i = 0; i < 2; ++i)
	{
		copyProbe(&k2jw.h2[i], jw->h2[i]);
		copyProbe(&k3jw.h2[i], jw->h2[i]);
		copyProbe(&k4jw.h2[i], jw->h2[i]);
	}

//	std::cout << "k2" << std::endl;
//	std::cout << "Printing position of jw" << std::endl;
//	std::cout << jw->pos.x << " " << jw->pos.y << std::endl;
//	std::cout << "Printing velocity of jw" << std::endl;
//	std::cout << jw->vel.x << " " << jw->vel.y << std::endl;
//	std::cout << "Printing position of v1" << std::endl;
//	std::cout << k2jw.h2[0].pos.x << " " << k2jw.h2[0].pos.y << std::endl;
//	std::cout << "Printing velocity of v1" << std::endl;
//	std::cout << k2jw.h2[0].vel.x << " " << k2jw.h2[0].vel.y << std::endl;
//	std::cout << "Printing position of v2" << std::endl;
//	std::cout << k2jw.h2[1].pos.x << " " << k2jw.h2[1].pos.y << std::endl;
//	std::cout << "Printing velocity of v2" << std::endl;
//	std::cout << k2jw.h2[1].vel.x << " " << k2jw.h2[1].vel.y << std::endl << std::endl;
//
//	std::cout << "k3" << std::endl;
//	std::cout << "Printing position of jw" << std::endl;
//	std::cout << jw->pos.x << " " << jw->pos.y << std::endl;
//	std::cout << "Printing velocity of jw" << std::endl;
//	std::cout << jw->vel.x << " " << jw->vel.y << std::endl;
//	std::cout << "Printing position of v1" << std::endl;
//	std::cout << "Printing position of v1" << std::endl;
//	std::cout << k3jw.h2[0].pos.x << " " << k3jw.h2[0].pos.y << std::endl;
//	std::cout << "Printing velocity of v1" << std::endl;
//	std::cout << k3jw.h2[0].vel.x << " " << k3jw.h2[0].vel.y << std::endl;
//	std::cout << "Printing position of v2" << std::endl;
//	std::cout << k3jw.h2[1].pos.x << " " << k3jw.h2[1].pos.y << std::endl;
//	std::cout << "Printing velocity of v2" << std::endl;
//	std::cout << k3jw.h2[1].vel.x << " " << k3jw.h2[1].vel.y << std::endl << std::endl;
//
//	std::cout << "k4" << std::endl;
//	std::cout << "Printing position of jw" << std::endl;
//	std::cout << jw->pos.x << " " << jw->pos.y << std::endl;
//	std::cout << "Printing velocity of jw" << std::endl;
//	std::cout << jw->vel.x << " " << jw->vel.y << std::endl;
//	std::cout << "Printing position of v1" << std::endl;
//	std::cout << k4jw.h2[0].pos.x << " " << k4jw.h2[0].pos.y << std::endl;
//	std::cout << "Printing velocity of v1" << std::endl;
//	std::cout << k4jw.h2[0].vel.x << " " << k4jw.h2[0].vel.y << std::endl;
//	std::cout << "Printing position of v2" << std::endl;
//	std::cout << k4jw.h2[1].pos.x << " " << k4jw.h2[1].pos.y << std::endl;
//	std::cout << "Printing velocity of v2" << std::endl;
//	std::cout << k4jw.h2[1].vel.x << " " << k4jw.h2[1].vel.y << std::endl << std::endl;

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K3 AFTER STEP" << std::endl;
//	std::cout << "Printing k3 object positions" << std::endl;
//	std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
//	std::cout << "k3Earth: v: " << k3Earth.vel.x << " " << k3Earth.vel.y << std::endl;
//	std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;
//	std::cout << "k3JamesWebb: v: " << k3jw.vel.x << " " << k3jw.vel.y << std::endl << std::endl << std::endl << std::endl;

	// k1 step

	// Calculate object accelerations and store in k1Accel
	calculateAccels(earth, sun, jw, NUMPERTURBED, &k1EarthAccel, &k1jwAccel, k1PerturbedAccels, k1h2Accels);

	// Unit testing
	// std::cout << "K1 STEP" << std::endl;
	// calcAccelTest(earth, sun, jw, NUMPERTURBED, &k1EarthAccel, &k1jwAccel, k1PerturbedAccels, k1h2Accels);

	// Extract velocities and store in k1
	extractVel(earth, sun, jw, NUMPERTURBED, &k1EarthVel, &k1jwVel, k1PerturbedVels, k1h2Vels);
	// extractVelTest(earth, sun, jw, NUMPERTURBED, &k1EarthVel, &k1jwVel, k1PerturbedVels, k1h2Vels);

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K3 AFTER STEP" << std::endl;
//	std::cout << "Printing k3 object positions" << std::endl;
//	std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
//	std::cout << "k3Earth: v: " << k3Earth.vel.x << " " << k3Earth.vel.y << std::endl;
//	std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;
//	std::cout << "k3JamesWebb: v: " << k3jw.vel.x << " " << k3jw.vel.y << std::endl << std::endl << std::endl << std::endl;

//	 DEBUGGING
//	std::cout << "Printing k1Accels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k1Accel[i].x << " y: " << k1Accel[i].y << std::endl;
//	}
//	std::cout << std::endl << "Printing k1Vels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k1Vel[i].x << " y: " << k1Vel[i].y << std::endl;
//	}

	// k2 Step

	// Move k2 object copies dt / 2 using k1.
//	std::cout << std::endl << "DEBUGGING K2 BEFORE STEP" << std::endl;
//	std::cout << "Printing k2 object positions" << std::endl;
//	std::cout << "k2Earth: x: " << k2Earth.pos.x << " " << k2Earth.pos.y << std::endl;
//	std::cout << "k2JamesWebb: x: " << k2jw.pos.x << " " << k2jw.pos.y << std::endl;
//	std::cout << std::endl << "k2" << std::endl;

	// std::cout << "K2 STEP" << std::endl;
	rk4Step(&k2Earth, &k2jw, NUMPERTURBED, &k1EarthAccel, &k1jwAccel, k1PerturbedAccels, k1h2Accels, dt / 2.0);
	// rk4StepTest(&k2Earth, &k2jw, NUMPERTURBED, &k1EarthAccel, &k1jwAccel, k1PerturbedAccels, k1h2Accels, dt / 2.0);

	// Calculate and store accelerations in k2
	calculateAccels(&k2Earth, sun, &k2jw, NUMPERTURBED, &k2EarthAccel, &k2jwAccel, k2PerturbedAccels, k2h2Accels);
	// calcAccelTest(&k2Earth, sun, &k2jw, NUMPERTURBED, &k2EarthAccel, &k2jwAccel, k2PerturbedAccels, k2h2Accels);


	// Extract velocities and store in k2
	extractVel(&k2Earth, sun, &k2jw, NUMPERTURBED, &k2EarthVel, &k2jwVel, k2PerturbedVels, k2h2Vels);
	// extractVelTest(&k2Earth, sun, &k2jw, NUMPERTURBED, &k2EarthVel, &k2jwVel, k2PerturbedVels, k2h2Vels);
	// DEBUGGING
//	std::cout << "Printing k2Accels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k2Accel[i].x << " y: " << k2Accel[i].y << std::endl;
//	}
//	std::cout << std::endl << "Printing k2Vels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k2Vel[i].x << " y: " << k2Vel[i].y << std::endl;
//	}
//
//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K3 AFTER STEP" << std::endl;
//	std::cout << "Printing k3 object positions" << std::endl;
//	std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
//	std::cout << "k3Earth: v: " << k3Earth.vel.x << " " << k3Earth.vel.y << std::endl;
//	std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;
//	std::cout << "k3JamesWebb: v: " << k3jw.vel.x << " " << k3jw.vel.y << std::endl << std::endl << std::endl << std::endl;
//
//	std::cout << std::endl << std::endl << std::endl << std::endl << "DEBUGGING K3 BEFORE STEP" << std::endl;
//	std::cout << "Printing k3 object positions" << std::endl;
//	std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
//	std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;

	// k3 Step
	// std::cout << "K3 STEP" << std::endl;

	// Move k3 object copies dt / 2 using k2
//	std::cout << std::endl << "k3" << std::endl;
	rk4Step(&k3Earth, &k3jw, NUMPERTURBED, &k2EarthAccel, &k2jwAccel, k2PerturbedAccels, k2h2Accels, dt / 2.0);
	// rk4StepTest(&k3Earth, &k3jw, NUMPERTURBED, &k2EarthAccel, &k2jwAccel, k2PerturbedAccels, k2h2Accels, dt / 2.0);

	// Calculate and store accelerations in k3
	calculateAccels(&k3Earth, sun, &k3jw, NUMPERTURBED, &k3EarthAccel, &k3jwAccel, k3PerturbedAccels, k3h2Accels);
	// calcAccelTest(&k3Earth, sun, &k3jw, NUMPERTURBED, &k3EarthAccel, &k3jwAccel, k3PerturbedAccels, k3h2Accels);

	// Extract velocities and store in k3
	extractVel(&k3Earth, sun, &k3jw, NUMPERTURBED, &k3EarthVel, &k3jwVel, k3PerturbedVels, k3h2Vels);
	// extractVelTest(&k3Earth, sun, &k3jw, NUMPERTURBED, &k3EarthVel, &k3jwVel, k3PerturbedVels, k3h2Vels);

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K3 AFTER STEP" << std::endl;
//	std::cout << "Printing k3 object positions" << std::endl;
//	std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
//	std::cout << "k3Earth: v: " << k3Earth.vel.x << " " << k3Earth.vel.y << std::endl;
//	std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;
//	std::cout << "k3JamesWebb: v: " << k3jw.vel.x << " " << k3jw.vel.y << std::endl << std::endl << std::endl << std::endl;

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;

	// DEBUGGING
//	std::cout << "Printing k3Accels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
////		std::cout << std::endl << "DEBUGGING K3 DURING LOOP EXTRACT" << std::endl;
////		std::cout << "Printing k3 object positions" << std::endl;
////		std::cout << "k3Earth: x: " << k3Earth.pos.x << " " << k3Earth.pos.y << std::endl;
////		std::cout << "k3Earth: v: " << k3Earth.vel.x << " " << k3Earth.vel.y << std::endl;
////		std::cout << "k3JamesWebb: x: " << k3jw.pos.x << " " << k3jw.pos.y << std::endl;
////		std::cout << "k3JamesWebb: v: " << k3jw.vel.x << " " << k3jw.vel.y << std::endl;
//		std::cout << "i: " << i << " x:" << k3Accel[i].x << " y: " << k3Accel[i].y << std::endl;
//	}
//	std::cout << std::endl << "Printing k3Vels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k3Vel[i].x << " y: " << k3Vel[i].y << std::endl;
//	}

	// k4 Step
	// std::cout << "K4 STEP" << std::endl;

	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
//	std::cout << std::endl << "k4" << std::endl;
	rk4Step(&k4Earth, &k4jw, NUMPERTURBED, &k3EarthAccel, &k3jwAccel, k3PerturbedAccels, k3h2Accels, dt);
	// rk4StepTest(&k4Earth, &k4jw, NUMPERTURBED, &k3EarthAccel, &k3jwAccel, k3PerturbedAccels, k3h2Accels, dt);

	// Calculate and store accelerations in k4
	calculateAccels(&k4Earth, sun, &k4jw, NUMPERTURBED, &k4EarthAccel, &k4jwAccel, k4PerturbedAccels, k4h2Accels);
	// calcAccelTest(&k4Earth, sun, &k4jw, NUMPERTURBED, &k4EarthAccel, &k4jwAccel, k4PerturbedAccels, k4h2Accels);

	// Extract velocities and store in k4
	extractVel(&k4Earth, sun, &k4jw, NUMPERTURBED, &k4EarthVel, &k4jwVel, k4PerturbedVels, k4h2Vels);
	// extractVelTest(&k4Earth, sun, &k4jw, NUMPERTURBED, &k4EarthVel, &k4jwVel, k4PerturbedVels, k4h2Vels);

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;


//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;


//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;


	// DEBUGGING
//	std::cout << "Printing k4Accels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k4Accel[i].x << " y: " << k4Accel[i].y << std::endl;
//	}
//	std::cout << std::endl << "Printing k4Vels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << k4Vel[i].x << " y: " << k4Vel[i].y << std::endl;
//	}

	//rk4 step
//	std::cout << "FINAL STEP" << std::endl << std::endl;

	// Fill up weighted k's

	// Earth:
	// Accelerations
	weightedEarthAccel.x = k1EarthAccel.x + 2 * k2EarthAccel.x + 2 * k3EarthAccel.x + k4EarthAccel.x;
	weightedEarthAccel.y = k1EarthAccel.y + 2 * k2EarthAccel.y + 2 * k3EarthAccel.y + k4EarthAccel.y;

	// Velocities
	weightedEarthVel.x = k1EarthVel.x + 2 * k2EarthVel.x + 2 * k3EarthVel.x + k4EarthVel.x;
	weightedEarthVel.y = k1EarthVel.y + 2 * k2EarthVel.y + 2 * k3EarthVel.y + k4EarthVel.y;

//	std::cout << "Weighted velocity: x: " << weightedEarthVel.x << " y: " << weightedEarthVel.y << std::endl;

	// James Webb:
	// Accelerations
	weightedjwAccel.x = k1jwAccel.x + 2 * k2jwAccel.x + 2 * k3jwAccel.x + k4jwAccel.x;
	weightedjwAccel.y = k1jwAccel.y + 2 * k2jwAccel.y + 2 * k3jwAccel.y + k4jwAccel.y;


//	std::cout << "James Webb" << std::endl;
//	std::cout << "Weighted acceleration: x: " << weightedjwAccel.x << " y: " << weightedjwAccel.y << std::endl;

	// Velocities
	weightedjwVel.x = k1jwVel.x + 2 * k2jwVel.x + 2 * k3jwVel.x + k4jwVel.x;
	weightedjwVel.y = k1jwVel.y + 2 * k2jwVel.y + 2 * k3jwVel.y + k4jwVel.y;

//	std::cout << "Weighted Velocity: x: " << weightedjwVel.x << " y: " << weightedjwVel.y << std::endl;

	// Perturbations:
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Accelerations
		weightedPerturbedAccels[i].x = k1PerturbedAccels[i].x + 2 * k2PerturbedAccels[i].x + 2 * k3PerturbedAccels[i].x + k4PerturbedAccels[i].x;
		weightedPerturbedAccels[i].y = k1PerturbedAccels[i].y + 2 * k2PerturbedAccels[i].y + 2 * k3PerturbedAccels[i].y + k4PerturbedAccels[i].y;

//		std::cout << "Perturbed Probe: " << i << std::endl;
//		std::cout << "Weighted acceleration: x: " << weightedPerturbedAccels[i].x << " y: " << weightedPerturbedAccels[i].y << std::endl;

		// Velocities
		weightedPerturbedVels[i].x = k1PerturbedVels[i].x + 2 * k2PerturbedVels[i].x + 2 * k3PerturbedVels[i].x + k4PerturbedVels[i].x;
		weightedPerturbedVels[i].y = k1PerturbedVels[i].y + 2 * k2PerturbedVels[i].y + 2 * k3PerturbedVels[i].y + k4PerturbedVels[i].y;

//		std::cout << "Weighted velocity: x: " << weightedPerturbedVels[i].x << " y: " << weightedPerturbedVels[i].y << std::endl;
	}

	// h2:
	for (int i = 0; i < 2; ++i)
	{
		// Accelerations
		weightedh2Accels[i].x = k1h2Accels[i].x + 2 * k2h2Accels[i].x + 2 * k3h2Accels[i].x + k4h2Accels[i].x;
		weightedh2Accels[i].y = k1h2Accels[i].y + 2 * k2h2Accels[i].y + 2 * k3h2Accels[i].y + k4h2Accels[i].y;

//		std::cout << "h2 Probe: " << i << std::endl;
//		std::cout << "Weighted acceleration: x: " << weightedh2Accels[i].x << " y: " << weightedh2Accels[i].y << std::endl;

		// Velocities
		weightedh2Vels[i].x = k1h2Vels[i].x + 2 * k2h2Vels[i].x + 2 * k3h2Vels[i].x + k4h2Vels[i].x;
		weightedh2Vels[i].y = k1h2Vels[i].y + 2 * k2h2Vels[i].y + 2 * k3h2Vels[i].y + k4h2Vels[i].y;

//		std::cout << "Weighted velocity: x: " << weightedh2Vels[i].x << " y: " << weightedh2Vels[i].y << std::endl;
	}
//	std::cout << std::endl << std::endl;

	// Step real objects forward dt / 6 using weighted k's
	// Earth:
	// Velocity

//	std::cout << "Earth" << std::endl;
//	std::cout << "Old velocities: vx: " << earth->vel.x << " vy: " << earth->vel.y << std::endl;

	earth->vel.x += weightedEarthAccel.x * dt / 6.0;
	earth->vel.y += weightedEarthAccel.y * dt / 6.0;

//	std::cout << "New velocities: vx: " << earth->vel.x << " vy: " << earth->vel.y << std::endl;
//	std::cout << "Old positions: x: " << earth->pos.x << " y: " << earth->pos.y << std::endl;

	// Position
	earth->pos.x += weightedEarthVel.x * dt / 6.0;
	earth->pos.y += weightedEarthVel.y * dt / 6.0;

//	std::cout << "New positions: x: " << earth->pos.x << " y: " << earth->pos.y << std::endl;

	// James Webb:
	// Velocity

//	std::cout << "James Webb" << std::endl;
//	std::cout << "Old velocities: vx: " << jw->vel.x << " vy: " << jw->vel.y << std::endl;

	jw->vel.x += weightedjwAccel.x * dt / 6.0;
	jw->vel.y += weightedjwAccel.y * dt / 6.0;

//	std::cout << "New velocities: vx: " << jw->vel.x << " vy: " << jw->vel.y << std::endl;
//	std::cout << "Old positions: x: " << jw->pos.x << " y: " << jw->pos.y << std::endl;

	// Position:
	jw->pos.x += weightedjwVel.x * dt / 6.0;
	jw->pos.y += weightedjwVel.y * dt / 6.0;

//	std::cout << "New positions: x: " << jw->pos.x << " y: " << jw->pos.y << std::endl;

	// Perturbed Probes:
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Velocity
//		std::cout << "Perturbed probe: " << i << std::endl;
//		std::cout << "Old velocities: vx: " << jw->perturbed[i].vel.x << " vy: " << jw->perturbed[i].vel.y << std::endl;

		jw->perturbed[i].vel.x += weightedPerturbedAccels[i].x * dt / 6.0;
		jw->perturbed[i].vel.y += weightedPerturbedAccels[i].y * dt / 6.0;

//		std::cout << "New velocities: vx: " << jw->perturbed[i].vel.x << " vy: " << jw->perturbed[i].vel.y << std::endl;
//		std::cout << "Old positions: x: " << jw->perturbed[i].pos.x << " y: " << jw->perturbed[i].pos.y << std::endl;

		// Position
		jw->perturbed[i].pos.x += weightedPerturbedVels[i].x * dt / 6.0;
		jw->perturbed[i].pos.y += weightedPerturbedVels[i].y * dt / 6.0;

//		std::cout << "New positions: x: " << jw->perturbed[i].pos.x << " y: " << jw->perturbed[i].pos.y << std::endl;
	}

	// h2:
	for (int i = 0; i < 2; ++i)
	{
		// Velocity

//		std::cout << "h2 probe: " << i << std::endl;
//		std::cout << "Old velocities: vx: " << jw->h2[i].vel.x << " vy: " << jw->h2[i].vel.y << std::endl;

		jw->h2[i].vel.x += weightedh2Accels[i].x * dt / 6.0;
		jw->h2[i].vel.y += weightedh2Accels[i].y * dt / 6.0;

//		std::cout << "New velocities: vx: " << jw->h2[i].vel.x << " vy: " << jw->h2[i].vel.y << std::endl;
//		std::cout << "Old positions: x: " << jw->h2[i].pos.x << " y: " << jw->h2[i].pos.y << std::endl;

		// Positions
		jw->h2[i].pos.x += weightedh2Vels[i].x * dt / 6.0;
		jw->h2[i].pos.y += weightedh2Vels[i].y * dt / 6.0;

//		std::cout << "New positions: x: " << jw->h2[i].pos.x << " y: " << jw->h2[i].pos.y << std::endl;
	}
//	std::cout << std::endl << std::endl << "--------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

	// rk4Step(earth, jw, NUMPERTURBED, &weightedEarthAccel, &weightedjwAccel, weightedPerturbedAccels, weightedh2Accels, dt / 6.0);
	//rk4StepTest(earth, jw, NUMPERTURBED, &weightedEarthAccel, &weightedjwAccel, weightedPerturbedAccels, weightedh2Accels, dt / 6.0);

	// DEBUGGING
//	std::cout << "Printing weightedKAccels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << weightedKAccel[i].x << " y: " << weightedKAccel[i].y << std::endl;
//	}
//	std::cout << std::endl << "Printing weightedKVels" << std::endl;
//	for (int i = 0; i < NUMPERTURBED + 4; ++i)
//	{
//		std::cout << "i: " << i << " x:" << weightedKVel[i].x << " y: " << weightedKVel[i].y << std::endl;
//	}

	// Free stuff

	// Free perturbed probes
	free(k2jw.perturbed);
	free(k3jw.perturbed);
	free(k4jw.perturbed);

	// Free h2 probes
	free(k2jw.h2);
	free(k3jw.h2);
	free(k4jw.h2);

	// Free acceleration and velocity arrays

	// Perturbed Probes:
	// Accelerations
	free(k1PerturbedAccels);
	free(k2PerturbedAccels);
	free(k3PerturbedAccels);
	free(k4PerturbedAccels);
	free(weightedPerturbedAccels);

	// Velocities
	free(k1PerturbedVels);
	free(k2PerturbedVels);
	free(k3PerturbedVels);
	free(k4PerturbedVels);
	free(weightedPerturbedVels);

	// h2:
	// Accelerations
	free(k1h2Accels);
	free(k2h2Accels);
	free(k3h2Accels);
	free(k4h2Accels);
	free(weightedh2Accels);

	// Velocities
	free(k1h2Vels);
	free(k2h2Vels);
	free(k3h2Vels);
	free(k4h2Vels);
	free(weightedh2Vels);
}


// Moves object Copies forward using kAccel and kVel
void singleProbeRK4Step(Probe *jwCopy, Planet *earthCopy, int NUMPERTURBED, vector2 *accels, vector2 *vels, double dt)
{
	// Earth
	earthCopy->vel.x += accels[0].x * dt; earthCopy->vel.y += accels[0].y * dt;
	earthCopy->pos.x += vels[0].x * dt; earthCopy->pos.y += vels[0].y * dt;

	// James webb probe
	jwCopy->vel.x += accels[1].x * dt; jwCopy->vel.y += accels[1].y * dt;
	jwCopy->pos.x += vels[1].x * dt; jwCopy->pos.y += vels[1].y * dt;

	// Perturbations
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		jwCopy->perturbed[i].vel.x += accels[i + 2].x * dt; jwCopy->perturbed[i].vel.y += accels[i + 2].y * dt;
		jwCopy->perturbed[i].pos.x += vels[i + 2].x * dt; jwCopy->perturbed[i].pos.y += vels[i + 2].y * dt;
	}
//
//	// h2 Probes
	for (int i = 0; i < 2; ++i)
	{
		jwCopy->h2[i].vel.x += accels[i + 2 + NUMPERTURBED].x * dt; jwCopy->h2[i].vel.y += accels[i + 2 + NUMPERTURBED].y * dt;
		jwCopy->h2[i].pos.x += vels[i + 2 + NUMPERTURBED].x * dt; jwCopy->h2[i].pos.y += vels[i + 2 + NUMPERTURBED].y * dt;
	}
}

// Moves Objects forward using kAccel and velocity
void rk4Step(Planet *earth, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels, double dt)
{
	// Push velocities forward using k Accel, and then push positions forward using current velocity

	// Earth:
	// Velocities
	earth->vel.x += earthAccel->x * dt;
	earth->vel.y += earthAccel->y * dt;

	// Positions
	earth->pos.x += earth->vel.x * dt;
	earth->pos.y += earth->vel.y * dt;

	// James Webb:
	// Velocities
	jw->vel.x += jwAccel->x * dt;
	jw->vel.y += jwAccel->y * dt;

	// Positions
	jw->pos.x += jw->vel.x * dt;
	jw->pos.y += jw->vel.y * dt;

	// Perturbations:
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Velocities
		jw->perturbed[i].vel.x += perturbedAccels[i].x * dt;
		jw->perturbed[i].vel.y += perturbedAccels[i].y * dt;

		// Positions
		jw->perturbed[i].pos.x += jw->perturbed[i].vel.x * dt;
		jw->perturbed[i].pos.y += jw->perturbed[i].vel.y * dt;
	}

	// h2:
	for (int i = 0; i < 2; ++i)
	{
		// Velocities
		jw->h2[i].vel.x += h2Accels[i].x * dt;
		jw->h2[i].vel.y += h2Accels[i].y * dt;

		// Positions
		jw->h2[i].pos.x += jw->h2[i].vel.x * dt;
		jw->h2[i].pos.y += jw->h2[i].vel.y * dt;
	}
}

// Stores velocities in kVels
void extractVel(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthVel, vector2 *jwVel, vector2 *perturbedVels, vector2 *h2Vels)
{
	// Earth
	earthVel->x = earth->vel.x;
	earthVel->y = earth->vel.y;

	// James Webb
	jwVel->x = jw->vel.x;
	jwVel->y = jw->vel.y;

	// Perturbed
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		perturbedVels[i].x = jw->perturbed[i].vel.x;
		perturbedVels[i].y = jw->perturbed[i].vel.y;
	}

	// h2
	for (int i = 0; i < 2; ++i)
	{
		h2Vels[i].x = jw->h2[i].vel.x;
		h2Vels[i].y = jw->h2[i].vel.y;
	}
}

// Stores velocities in kVel
void extractkVelSingleProbe(Planet earth, Probe jw, int NUMPERTURBED, vector2 *vels)
{
//	std::cout << std::endl << std::endl << "DEBUGGING EXTRACT VELOCITY" << std::endl;
//	std::cout << "Printing old earth velocity: " << earth.vel.x << " " << earth.vel.y << std::endl;

	vels[0].x = earth.vel.x; vels[0].y = earth.vel.y;
//	std::cout << "Printing new earth velocity: " << earth.vel.x << " " << earth.vel.y << std::endl;
//	std::cout << "Printing old jw velocity: " << jw.vel.x << " " << jw.vel.y << std::endl;
	vels[1].x = jw.vel.x; vels[1].y = jw.vel.y;

//	std::cout << "Printing new jw velocity: " << jw.vel.x << " " << jw.vel.y << std::endl;

	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		vels[i + 2].x = jw.perturbed[i].vel.x; vels[i + 2].y = jw.perturbed[i].vel.y;
	}
	for (int i = 0; i < 2; ++i)
	{
		vels[i + 2 + NUMPERTURBED].x = jw.h2[i].vel.x; vels[i + 2 + NUMPERTURBED].y = jw.h2[i].vel.y;
	}
}

// Calculates and stores gravitational accelerations on earth, james webb, perturbations, and h2
void calculateAccels(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels)
{
	double G = 6.6743e-11;

	// Sun on earth:
	vector2 r = {earth->pos.x - sun->pos.x, earth->pos.y - sun->pos.y}; 	// Vector pointing from sun to earth
	double magr = mag(r);														// Magnitude of that vector
	vector2 rUnit = unitVector(r);												// Direction of the vector
	double g = (-1 * G * sun->m) / (magr * magr);								// Newtons law of gravity to get force magnitude
	// Multiply by unit vector for direction
	earthAccel->x = g * rUnit.x;
	earthAccel->y = g * rUnit.y;

	// James Webb Telescope: Add up accelerations from Sun and Earth: Same process

	// Sun on James Webb
	r = {jw->pos.x - sun->pos.x, jw->pos.y - sun->pos.y};
	magr = mag(r);
	rUnit = unitVector(r);
	g = (-1 * G * sun->m) / (magr * magr);
	jwAccel->x = g * rUnit.x;
	jwAccel->y = g * rUnit.y;

	// Earth on James Webb
	r = {jw->pos.x - earth->pos.x, jw->pos.y - earth->pos.y};
	magr = mag(r);
	rUnit = unitVector(r);
	g = (-1 * G * earth->m) / (magr * magr);
	// Add to accelerations from the Sun
	jwAccel->x += g * rUnit.x;
	jwAccel->y += g * rUnit.y;

	// Perturbations
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Sun on perturbation
		r = {jw->perturbed[i].pos.x - sun->pos.x, jw->perturbed[i].pos.y - sun->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * sun->m) / (magr * magr);
		perturbedAccels[i].x = g * rUnit.x;
		perturbedAccels[i].y = g * rUnit.y;

		// Earth on perturbation
		r = {jw->perturbed[i].pos.x - earth->pos.x, jw->perturbed[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		perturbedAccels[i].x += g * rUnit.x;
		perturbedAccels[i].y += g * rUnit.y;
	}

	// h2 Probes
	for (int i = 0; i < 2; ++i)
	{
		// Sun on h2 probe
		r = {jw->h2[i].pos.x - sun->pos.x, jw->h2[i].pos.y - sun->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * sun->m) / (magr * magr);
		h2Accels[i].x = g * rUnit.x;
		h2Accels[i].y = g * rUnit.y;

		// Earth on h2 probe
		r = {jw->h2[i].pos.x - earth->pos.x, jw->h2[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		h2Accels[i].x += g * rUnit.x;
		h2Accels[i].y += g * rUnit.y;
	}
}

void extractkAccelSingleProbe(Planet earth, Planet sun, Probe jw, int NUMPERTURBED, vector2 *accels)
{
	double G = 6.6743e-11;

	// Sun: 0, Earth: 1, James Webb: 2, Perturbations: 3, h2 Probes: 4
	// Sun on earth: Sun on Earth
	vector2 r01 = {earth.pos.x - sun.pos.x, earth.pos.y - sun.pos.y};
	double magr01 = mag(r01);
	vector2 r01Unit = unitVector(r01);
	double g01 = (-1 * G * sun.m) / (magr01 * magr01);
	accels[0].x = g01 * r01Unit.x; accels[0].y = g01 * r01Unit.y;

//	std::cout << "Sun on earth acceleration is: " << accels[0].x << " " << accels[0].y << std::endl;
//	std::cout << "Earth: mass: " << earth.m << " Sun: mass: " << sun.m << std::endl;
//	std::cout << "Dist: " << magr01 << std::endl;
//	std::cout << std::endl << "Acceleration Sun on Earth" << std::endl;
//	std::cout << "Sun on Earth using doubles: Force: " << g01 << " Acceleration: x: " << accels[0].x << " y: " << accels[0].y << std::endl;

	// Sun on james webb telescope

//	std::cout << "BEFORE CALCULATING SUN ON JAMES WEBB" << std::endl;
//	std::cout << "Printing James webb data" << std::endl;
//	std::cout << "k2JamesWebb: x: " << jw.pos.x << " " << jw.pos.y << std::endl;

	// Sun on james webb telescope
	vector2 r02 = {jw.pos.x - sun.pos.x, jw.pos.y - sun.pos.y};
	double magr02 = mag(r02);
//	std::cout << "Printing james webb dist from sun: " << magr02 << std::endl;
	vector2 r02Unit = unitVector(r02);
//	std::cout << "Printing unit vector: x: " << r02Unit.x << " y: " << r02Unit.y << std::endl;
	double g02 = (-1 * G * sun.m) / (magr02 * magr02);
//	std::cout << "Printing g02: x: " << g02 << std::endl;
	vector2 res02 = {g02 * r02Unit.x, g02 * r02Unit.y};

//	std::cout << "Sun on james webb telescope is: " << g02 << std::endl;

	// DEBUGGING Earth on James Webb
//	std::cout << "Printing sun accelerations on jw: x: " << accels[1].x << " " << accels[1].y << std::endl;
	vector2 r12 = {jw.pos.x - earth.pos.x, jw.pos.y - earth.pos.y};
	double magr12 = mag(r12);
	vector2 r12Unit = unitVector(r12);
//	std::cout << "Printing r12: x: " << r12.x << " " << r12.y << std::endl;
 	double g12 = (-1 * G * earth.m) / (magr12 * magr12);
	vector2 res12 = {g12 * r12Unit.x, g12 * r12Unit.y};
	accels[1].x = res02.x + res12.x; accels[1].y = res02.y + res12.y;

	// Perturbations
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Sun on perturbed
		vector2 r03 = {jw.perturbed[i].pos.x - sun.pos.x, jw.perturbed[i].pos.y - sun.pos.y};
		double magr03 = mag(r03);
		vector2 r03Unit = unitVector(r03);
		double g03 = (-1 * G * sun.m) / (magr03 * magr03);
		vector2 res03 = {g03 * r03Unit.x, g03 * r03Unit.y};

		// Earth on perturbed
		vector2 r13 = {jw.perturbed[i].pos.x - earth.pos.x, jw.perturbed[i].pos.y - earth.pos.y};
		double magr13 = mag(r13);
		vector2 r13Unit = unitVector(r13);
		double g13 = (-1 * G * earth.m) / (magr13 * magr13);
		vector2 res13 = {g13 * r13Unit.x, g13 * r13Unit.y};

		accels[i + 2].x = res03.x + res13.x; accels[i + 2].y = res03.y + res13.y;
	}

	// h2
	for (int i = 0; i < 2; ++i)
	{
		// Sun on h2
		vector2 r04 = {jw.h2[i].pos.x - sun.pos.x, jw.h2[i].pos.y - sun.pos.y};
		double magr04 = mag(r04);
		vector2 r04Unit = unitVector(r04);
		double g04 = (-1 * G * sun.m) / (magr04 * magr04);
		vector2 res04 = {g04 * r04Unit.x, g04 * r04Unit.y};

		// Earth on h2
		vector2 r14 = {jw.h2[i].pos.x - earth.pos.x, jw.h2[i].pos.y - earth.pos.y};
		double magr14 = mag(r14);
		vector2 r14Unit = unitVector(r14);
		double g14 = (-1 * G * earth.m) / (magr14 * magr14);
		vector2 res14 = {g14 * r14Unit.x, g14 * r14Unit.y};

		accels[i + 2 + NUMPERTURBED].x = res04.x + res14.x; accels[i + 2 + NUMPERTURBED].y = res04.y + res14.y;
	}

//	std::cout << "Earth on james webb telescope is: " << g12 << std::endl;
//	std::cout << "Net on james webb telescope is: " << accels[1].x << " " << accels[1].y << std::endl;

//	std::cout << std::endl << std::endl << "Forces on james webb" << std::endl;
//	std::cout << "Magr12: " << magr12 << std::endl;
//	std::cout << "magr12 ^2: " << magr12 * magr12 << std::endl;
//	std::cout << "Numerator: " << top << std::endl;
//	std::cout << "Mass of earth: " << earth.m << std::endl;
//	std::cout << "Printing earth accelerations on jw: x: " << res12.x << " " << res12.y << std::endl;
//	std::cout << "Sun on James Webb using net: Force: " << g12 << " Acceleration: x: " << accels[1].x << " y: " << accels[1].y << std::endl;
//	std::cout << "Sun on James Webb using vector class: Force: " << gravity02 << " Acceleration: x: " << res.x << " y: " << res.y << std::endl;
//	std::cout << "Earth on James Webb using vector class: Force: " << gravity12 << " Acceleration: x: " << res12.x << " y: " << res12.y << std::endl;
//	std::cout << "Earth on James Webb using vector class: Force: " << g12 << " Acceleration: x: " << res12d[0] << " y: " << res12d[1] << std::endl;
//	std::cout << "Net on James Webb using vector class: x: " << accels[1].x << " y: " << accels[1].y << std::endl;


//	std::cout << "Earth on James webb: x: " << r12.x << " y: " << r12.y << std::endl;
//	std::cout << "Net accel: " << accels[1].x << " " << accels[1].y << std::endl;
//	std::cout << "Net accels: " << accels[1].x << " " << accels[1].y << std::endl;
//	std::cout << "Printing calculated gravitational acceleration: " << g12 << " Printing live calculated: " << g12 << std::endl;
//	std::cout << "Printing the unit vector rUnit: x: " << r12Unit.x << " y: " << r12Unit.y << std::endl;
//	std::cout << "The acceleration vector that was added: x: " << g12 * r12Unit.x << " y: " << g12 * r12Unit.y << std::endl << std::endl;
}
// Initializes the perturbations around a single probe
void initializePerturbations(Probe *jw, int NUMPERTURBED)
{
	// Initialize all perturbation probes: Last five: 0,1 are for h2 lyapunov: 2,3,4 are for h3 lyapunov
	jw->perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		copyProbe(&jw->perturbed[i], *jw);
	}

	// Read from file
	// Desktop path
	FILE *f = fopen("D:\\2dSunEarthChaos\\unitCircle.csv", "r");
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);
	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);
	// Read from buffer into perturbed probes
	int pointCount = 0;
	char arr[30];
	int arrCount = 0;
	int i = 0;

	while (pointCount < NUMPERTURBED)
	{
		if (buffer[i] == '\n')
		{
			arr[arrCount] = '\0';
			jw->perturbed[pointCount].pos.y = std::atof(arr) + jw->pos.y;
			pointCount++;
			arrCount = 0;
		}
		else if (buffer[i] == ',')
		{
			arr[arrCount] = '\0';
			jw->perturbed[pointCount].pos.x = std::atof(arr) + jw->pos.x;
			arrCount = 0;
		}
		else
		{
			arr[arrCount] = buffer[i];
			arrCount++;
		}
		++i;
	}
	free(buffer);
}
