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

	// Allocate h2 probes and copy over info:
	jw->h2 = (Probe *)malloc(2 * sizeof(Probe)); for (int i = 0; i < 2; ++i) { copyProbe(&jw->h2[i], *jw); }

//	std::cout << "INITIAL" << std::endl;
//	std::cout << "Printing position of JW" << std::endl;
//	std::cout << jw->pos.x << " " << jw->pos.y << std::endl;
//	std::cout << "Printing velocity of JW" << std::endl;
//	std::cout << jw->vel.x << " " << jw->vel.y << std::endl;
//	std::cout << "Printing position of v1" << std::endl;
//	std::cout << jw->h2[0].pos.x << " " << jw->h2[0].pos.y << std::endl;
//	std::cout << "Printing velocity of v1" << std::endl;
//	std::cout << jw->h2[0].vel.x << " " << jw->h2[0].vel.y << std::endl;
//	std::cout << "Printing position of v2" << std::endl;
//	std::cout << jw->h2[1].pos.x << " " << jw->h2[1].pos.y << std::endl;
//	std::cout << "Printing velocity of v2" << std::endl;
//	std::cout << jw->h2[1].vel.x << " " << jw->h2[1].vel.y << std::endl << std::endl;

	// Initialize h2 orthonormal vectors:
	jw->h2[0].pos.x += 1.0; jw->h2[1].pos.y += 1.0;

	// Storage for h1 (max) and h2 (min) lyapunov exponents
	double *lyapunovs = (double *)malloc(2 * sizeof(double));

	// Sum of log of relative separation every timestep for each perturbation: ORBITAL SEPARATION METHOD FOR CALCULATING LYAPUNOVS
	// Max lyapunov will be extracted from this at the end
	double *LfSums = (double *)malloc(NUMPERTURBED * sizeof(double)); for (int i = 0; i < NUMPERTURBED; ++i) { LfSums[i] = 0; }

	// h2 lyapunov will be extracted from this at the end
	double areaSum = 0;
	double initialArea = RADIUS * RADIUS;

	// Store data: Store all data of probe and perturbations every ITERATIONS step: Useful for plotting.
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
	vector2 **data = (vector2 **)malloc(n * sizeof(vector2 *));
	// 0: Earth, 1: probe, 2 - NUMPERTURBED + 2: perturbations, NUMPERTURBED + 2 - NUMPERTURBED + 3: h2
	for (int i = 0; i < n; ++i) { data[i] = (vector2 *)malloc((NUMPERTURBED + 10) * sizeof(vector2)); }
	// Lyapunov real time calculation
	double *lt = (double *)malloc(2 * sizeof(double));

	// Loop
	int c = 0;
	int dc = 0;
	bool t = true;
	while (c < totalSteps)
	{
		// Store position data of Probe, perturbations, and h2 orthogonals every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			// Earth
			data[dc][0].x = earth->pos.x; data[dc][0].y = earth->pos.y;

//			// James Webb
			data[dc][1].x = jw->pos.x; data[dc][1].y = jw->pos.y;

//			// Perturbations
//			for (int i = 0; i < NUMPERTURBED; ++i) { data[dc][i + 2].x = jw->perturbed[i].pos.x; data[dc][i + 2].y = jw->perturbed[i].pos.y; }

//			// h2
//			data[dc][NUMPERTURBED + 2].x = jw->h2[0].pos.x; data[dc][NUMPERTURBED + 2].y = jw->h2[0].pos.y; data[dc][NUMPERTURBED + 3].x = jw->h2[1].pos.x; data[dc][NUMPERTURBED + 3].y = jw->h2[1].pos.y;

//			// h2 normalized and offset for plotting
//			vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y}; vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};
//			vector2 u1 = unitVector(v1); vector2 u2 = unitVector(v2);
//			data[dc][NUMPERTURBED + 4].x = u1.x; data[dc][NUMPERTURBED + 4].y = u1.y; data[dc][NUMPERTURBED + 5].x = u2.x; data[dc][NUMPERTURBED + 5].y = u2.y;

			// Lyapunov exponent at current time step
//			extractLyapunovSingleProbe(LfSums, areaSum, NUMPERTURBED, lt, dt * (double)c);
			// h1 stored in x, h2 stored in y
//			data[dc][NUMPERTURBED + 6].x = lt[0]; data[dc][NUMPERTURBED + 6].y = lt[1];

			// Store max distance and area
//			double max = LfSums[0];
//			int maxInd = 0;
//			for (int i = 1; i < NUMPERTURBED; ++i)
//			{
//				if (LfSums[i] > max) { max = LfSums[i]; maxInd = i; }
//			}
//			vector2 maxV = {jw->perturbed[maxInd].pos.x - jw->pos.x, jw->perturbed[maxInd].pos.y - jw->pos.y};
//			// Store magnitude in x, area in y
//			data[dc][NUMPERTURBED + 7].x = mag(maxV);
//			vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y}; vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};
//			data[dc][NUMPERTURBED + 7].y = (mag(v1) * mag(v2) * sin(angle(v1, v2)));

			// Distance of l2 probe from original position over time: stability check
			vector2 sunToEarth = {data[dc][0].x - 149597870700.0, data[dc][0].y};
			vector2 dir = unitVector(sunToEarth);
			vector2 l2 = {1.5e9 * dir.x, 1.5e9 * dir.y};
//			std::cout << "unitx: " << dir.x << " unity: " << dir.y << std::endl;
			data[dc][NUMPERTURBED + 8].x = dir.x, data[dc][NUMPERTURBED + 8].y = dir.y;
			data[dc][NUMPERTURBED + 9].x = l2.x; data[dc][NUMPERTURBED + 9].x = l2.y;
			vector2 d = {l2.x - data[dc][1].x, l2.y - data[dc][1].y};
			std::cout << "L2 is at: " << l2.x << " " << l2.y << std::endl;
			std::cout << "Probe is at: " << data[dc][1].x << " " << data[dc][1].y << std::endl;
			std::cout << "Earth is at: " << data[dc][0].x << " " << data[dc][0].y << std::endl;
			std::cout << "Dist from l2: " << mag(d) << std::endl << std::endl;

			dc++;
		}
//		if (t)
//		{
//			if (jw->vel.x > 35000 || jw->vel.y > 35000)
//			{
//				std::cout << "Broken at iteration: " << c << std::endl;
//				std::cout << "Velocity is: " << jw->vel.x << " " << jw->vel.y << std::endl;
//				t = false;
//			}
//		}
		// Numerical Integrator: RK4
		singleProbeRK4(earth, sun, jw, NUMPERTURBED, dt);

//		std::cout << "Printing position of JW" << std::endl;
//		std::cout << jw->pos.x << " " << jw->pos.y << std::endl;
//		std::cout << "Printing velocity of JW" << std::endl;
//		std::cout << jw->vel.x << " " << jw->vel.y << std::endl;
//		std::cout << "Printing position of v1" << std::endl;
//		std::cout << jw->h2[0].pos.x << " " << jw->h2[0].pos.y << std::endl;
//		std::cout << "Printing velocity of v1" << std::endl;
//		std::cout << jw->h2[0].vel.x << " " << jw->h2[0].vel.y << std::endl;
//		std::cout << "Printing position of v2" << std::endl;
//		std::cout << jw->h2[1].pos.x << " " << jw->h2[1].pos.y << std::endl;
//		std::cout << "Printing velocity of v2" << std::endl;
//		std::cout << jw->h2[1].vel.x << " " << jw->h2[1].vel.y << std::endl;

		// Function for LfSums, Gram Schmidt2, and Renormalization of everything
		lyapunovChaosStuffSingle(jw, NUMPERTURBED, LfSums, RADIUS, &areaSum, initialArea);
		c++;
	}

	// Store latest info
	// Earth
//	data[n - 1][0].x = earth->pos.x; data[n - 1][0].y = earth->pos.y;
//	// James webb
//	data[n - 1][1].x = jw->pos.x; data[n - 1][1].y = jw->pos.y;

	// Storing perturbation info
//	for (int i = 0; i < NUMPERTURBED; ++i) { data[n - 1][i + 2] = {jw->perturbed[i].pos.x, jw->perturbed[i].pos.y}; }
//	data[n - 1][NUMPERTURBED + 2] = {jw->h2[0].pos.x, jw->h2[0].pos.y}; data[n - 1][NUMPERTURBED + 3] = {jw->h2[1].pos.x, jw->h2[1].pos.y};
//	vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y}; vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};
//	vector2 u1 = unitVector(v1); vector2 u2 = unitVector(v2);
//	data[n - 1][NUMPERTURBED + 4].x = u1.x; data[n - 1][NUMPERTURBED + 4].y = u1.y; data[n - 1][NUMPERTURBED + 5].x = u2.x; data[n - 1][NUMPERTURBED + 5].y = u2.y;

	// Storing lyapunov
//	extractLyapunovSingleProbe(LfSums, areaSum, NUMPERTURBED, lt, dt * (double)totalSteps);
//	data[n - 1][NUMPERTURBED + 6].x = lt[0]; data[n - 1][NUMPERTURBED + 6].y = lt[1];\

	// Dist and area
//	double max = LfSums[0];
//	int maxInd = 0;
//	for (int i = 1; i < NUMPERTURBED; ++i)
//	{
//		if (LfSums[i] > max) { max = LfSums[i]; maxInd = i; }
//	}
//	vector2 maxV = {jw->perturbed[maxInd].pos.x - jw->pos.x, jw->perturbed[maxInd].pos.y - jw->pos.y};
//	// Store magnitude in x, area in y
//	data[n - 1][NUMPERTURBED + 7].x = mag(maxV);
//	vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y}; vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};
//	data[n - 1][NUMPERTURBED + 7].y = (mag(v1) * mag(v2) * sin(angle(v1, v2)));

	// Distance of l2 probe from original position over time: stability check
//	vector2 sunToEarth = {data[dc][0].x - 149597870700.0, data[dc][0].y};
//	vector2 dir = unitVector(sunToEarth);
//	vector2 l2 = {1.5e9 * dir.x, 1.5e9 * dir.y};
//	data[n - 1][NUMPERTURBED + 8].x = dir.x, data[n - 1][NUMPERTURBED + 8].y = dir.y;
//	data[n - 1][NUMPERTURBED + 9].x = l2.x; data[n - 1][NUMPERTURBED + 9].x = l2.y;

	// Calculate and store lyapunov exponents
	extractLyapunovSingleProbe(LfSums, areaSum, NUMPERTURBED, lyapunovs, dt * (double)totalSteps);

	// Write data to file
//	dataToFile(data, n, NUMPERTURBED);

	// Free stuff
	for (int i = 0; i < n; ++i) { free(data[i]); } free(data); free(LfSums); free(jw->perturbed); free(jw->h2);
	// If using realtime
	free(lt);

	return lyapunovs;
}

// Writes position data to csv file for all objects in scene for plotting
void dataToFile(vector2 **data, int dataLen, int NUMPERTURBED)
{
	// File paths and names
	std::string fpath = "D:\\2dSunEarthChaos\\singleProbeData\\";
//	std::string f1 = "earthAndprobe.csv";
//	std::string f2 = "perturbations.csv";
//	std::string f3 = "gramSchmidtVectors.csv";
//	std::string f4 = "gramSchmidtProbes.csv";
//	std::string f5 = "NormalizedGramSchmidt.csv";
	std::string f6 = "leOverTime.csv";
//	std::string f7 = "distArea.csv";
	std::string f8 = "distFromL2.csv";
//	std::ofstream file1; std::ofstream file2; std::ofstream file3; std::ofstream file4; std::ofstream file5;
	std::ofstream file6;
//	std::ofstream file7;
	std::ofstream file8;
//	file1.open(fpath + f1); file2.open(fpath + f2); file3.open(fpath + f3); file4.open(fpath + f4); file5.open(fpath + f5);
	file6.open(fpath + f6);
//	file7.open(fpath + f7);
	file8.open(fpath + f8);
	// File headers
//	file1 << "ex,ey,jwx,jwy" << std::endl;
//	file3 << "v1x,v1y,v2x,v2y" << std::endl;
//	file4 << "v1x,v1y,v2x,v2y" << std::endl;
//	file5 << "u1x,u1y,u2x,u2y" << std::endl;
	file6 << "h1,h2" << std::endl;
//	file7 << "dist,area" << std::endl;
	file8 << "x,y,mag" << std::endl;
//	for (int i = 0; i < NUMPERTURBED; ++i)
//	{
//		file2 << "x" + std::to_string(i) << "," << "y" + std::to_string(i) << ",";
//	}
//	file2 << std::endl;

//	std::cout << "DataLen: " << dataLen << std::endl;
	// TODO: ADD THE DISTANCE DATA TO THE FILE
	for (int i = 0; i < dataLen; ++i)
	{
		// Earth
//		file1 << std::to_string(data[i][0].x) << "," << std::to_string(data[i][0].y) << ",";
		// James Webb
//		file1 << std::to_string(data[i][1].x) << "," << std::to_string(data[i][1].y) << std::endl;
//		for (int j = 0; j < NUMPERTURBED - 1; ++j)
//		{
//			file2 << std::to_string(data[i][j + 2].x) << "," << std::to_string(data[i][j + 2].y) << ",";
//		}
		// Last perturbation
//		file2 << std::to_string(data[i][NUMPERTURBED + 1].x) << "," << std::to_string(data[i][NUMPERTURBED + 1].y) << std::endl;
		// Two Gram schmidt vectors: x1,y1,x2,y2,u1x,u1y,u2x,u2y. Offset by probe's current position for plotting direction, with normalized option
//		file3 << std::to_string(data[i][NUMPERTURBED + 2].x - data[i][1].x) << "," << std::to_string(data[i][NUMPERTURBED + 2].y - data[i][1].y) << "," << std::to_string(data[i][NUMPERTURBED + 3].x - data[i][1].x) << "," << std::to_string(data[i][NUMPERTURBED + 3].y - data[i][1].y) << std::endl;
//		file5 << std::to_string(data[i][NUMPERTURBED + 4].x) << "," << std::to_string(data[i][NUMPERTURBED + 4].y) << "," << std::to_string(data[i][NUMPERTURBED + 5].x) << "," << std::to_string(data[i][NUMPERTURBED + 5].y) << std::endl;
//
//		// Two gram schmidt probes: no offset
//		file4 << std::to_string(data[i][NUMPERTURBED + 2].x) << "," << std::to_string(data[i][NUMPERTURBED + 2].y) << "," << std::to_string(data[i][NUMPERTURBED + 3].x) << "," << std::to_string(data[i][NUMPERTURBED + 3].y) << std::endl;

		// Lyapunov outputs
//		std::stringstream s1;
//		std::stringstream s2;
//		s1.precision(10); s2.precision(10);
//		s1 << data[i][NUMPERTURBED + 6].x; s2 << data[i][NUMPERTURBED + 6].y;
////
//		std::string out = s1.str() + "," + s2.str() + "\n";
//////		std::cout << out;
//		file6 << out;
		// Distance and area
//		std::stringstream s1;
//		std::stringstream s2;
//		s1.precision(10); s2.precision(10);
//		s1 << data[i][NUMPERTURBED + 7].x; s2 << data[i][NUMPERTURBED + 7].y;
//
//		std::string out = s1.str() + "," + s2.str() + "\n";
//		file7 << out;

		vector2 v = {data[i][NUMPERTURBED + 8].x, data[i][NUMPERTURBED + 8].y};
		std::stringstream s1;
		std::stringstream s2;
		std::stringstream s3;
		s1.precision(10); s2.precision(10); s3.precision(10);
		s1 << v.x; s2 << v.y; s3 << mag(v);
		std::string out = s1.str() + "," + s2.str() + "," + s3.str() + "\n";
		file8 << out;
	}

//	file1.close(); file2.close(); file3.close(); file4.close(); file5.close();
	file6.close();
//	file7.close();
	file8.close();
}
// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovSingleProbe(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t)
{
	// Calculate and store h1: max lyapunov
	double max = LfSums[0];
	for (int i = 1; i < NUMPERTURBED; ++i) { if (LfSums[i] > max) { max = LfSums[i]; } }
	lyapunovs[0] = max / t;
	// Calculate and store h2
	lyapunovs[1] = (areaSum / t) - lyapunovs[0];
}

// Handles the timestep processes for calculating h1 and h2 lyapunov exponents
void lyapunovChaosStuffSingle(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea)
{
	// h1: Max lyapunov
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Calculate Lf of each perturbation, and renormalize the position vector
		vector2 v = {jw->perturbed[i].pos.x - jw->pos.x, jw->perturbed[i].pos.y - jw->pos.y};
		// RADIUS = 1.0
		double dist = mag(v);
		double logDist = log(dist);

		// DEBUGGING
		// if (dist > maxMag) { maxMag = dist; }

		LfSums[i] += logDist;
		vector2 uv = unitVector(v);
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
	double logArea = log(area);
	(*areaSum) += logArea;
//	std::cout << std::setprecision(10) << area << std::endl;
//	std::cout << "magv1: " << mag(v1) << " magv2: " << mag(v2) << std::endl;

	// Reorthonormalize h2 vectors
	gramSchmidtRenormalizationh2Single(jw, NUMPERTURBED, v1, v2);
}

// Gram schmidt to re-orthonormalize unit vectors
void gramSchmidtRenormalizationh2Single(Probe *jw, int NUMPERTURBED, vector2 v1, vector2 v2)
{
	// Calculate orthonormal vectors e1 and e2
	vector2 e1 = unitVector(v1);
	double v2Dote1 = dotProd(v2, e1);
	vector2 y2 = {v2.x - (v2Dote1) * e1.x, v2.y - (v2Dote1) * e1.y};
	vector2 e2 = unitVector(y2);

//	std::cout << "e1: " << e1.x << " " << e1.y << std::endl;
//	std::cout << "e2: " << e2.x << " " << e2.y << std::endl << std::endl;

	// Assign new vectors
	jw->h2[0].pos.x = jw->pos.x + e1.x; jw->h2[0].pos.y = jw->pos.y + e1.y;
	jw->h2[1].pos.x = jw->pos.x + e2.x; jw->h2[1].pos.y = jw->pos.y + e2.y;
}
// RK4 for sun earth system with single james webb with perturbations
void singleProbeRK4(Planet *earth, Planet sun, Probe *jw, int NUMPERTURBED, double dt)
{
	// Initialize earth Copies
	Planet k2Earth, k3Earth, k4Earth;
	// Initialize jw Copies
	Probe k2jw, k3jw, k4jw;

	// Initialize k's: Index 0 is for earth, Index 1 is for probe, last 2 are for h2 probes
	vector2 *k1Accel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k2Accel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k3Accel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k4Accel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *weightedKAccel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2));
	vector2 *k1Vel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k2Vel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k3Vel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *k4Vel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2)); vector2 *weightedKVel = (vector2 *)malloc((NUMPERTURBED + 4) * sizeof(vector2));

	// Copy data into copies: planets, probes, perturbations, and h2 for all k
	copyPlanet(&k2Earth, *earth); copyPlanet(&k3Earth, *earth); copyPlanet(&k4Earth, *earth);
	copyProbe(&k2jw, *jw); copyProbe(&k3jw, *jw); copyProbe(&k4jw, *jw);
	k2jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe)); k3jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe)); k4jw.perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
	k2jw.h2 = (Probe *)malloc(2 * sizeof(Probe)); k3jw.h2 = (Probe *)malloc(2 * sizeof(Probe)); k4jw.h2 = (Probe *)malloc(2 * sizeof(Probe));
	for (int i = 0; i < NUMPERTURBED; ++i) { copyProbe(&k2jw.perturbed[i], jw->perturbed[i]); copyProbe(&k3jw.perturbed[i], jw->perturbed[i]); copyProbe(&k4jw.perturbed[i], jw->perturbed[i]); }
	for (int i = 0; i < 2; ++i)
	{
		copyProbe(&k2jw.h2[i], jw->h2[i]); copyProbe(&k3jw.h2[i], jw->h2[i]); copyProbe(&k4jw.h2[i], jw->h2[i]);
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
//	std::cout << std::endl << "k1" << std::endl;
	extractkAccelSingleProbe(*earth, sun, *jw, NUMPERTURBED, k1Accel);
	// Extract velocities from objects and store in k1Vel
	extractkVelSingleProbe(*earth, *jw, NUMPERTURBED, k1Vel);

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
	// Move k2 object copies dt / 2 using k1Accel and k1Vel.
//	std::cout << std::endl << "DEBUGGING K2 BEFORE STEP" << std::endl;
//	std::cout << "Printing k2 object positions" << std::endl;
//	std::cout << "k2Earth: x: " << k2Earth.pos.x << " " << k2Earth.pos.y << std::endl;
//	std::cout << "k2JamesWebb: x: " << k2jw.pos.x << " " << k2jw.pos.y << std::endl;
//	std::cout << std::endl << "k2" << std::endl;
	singleProbeRK4Step(&k2jw, &k2Earth, NUMPERTURBED, k1Accel, k1Vel, dt / 2.0);
	extractkAccelSingleProbe(k2Earth, sun, k2jw, NUMPERTURBED, k2Accel);
	extractkVelSingleProbe(k2Earth, k2jw, NUMPERTURBED, k2Vel);

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
	// Move k3 object copies dt / 2 using k2 Accelerations and k2 Velocities
//	std::cout << std::endl << "k3" << std::endl;
	singleProbeRK4Step(&k3jw, &k3Earth, NUMPERTURBED, k2Accel, k2Vel, dt / 2.0);
	extractkAccelSingleProbe(k3Earth, sun, k3jw, NUMPERTURBED, k3Accel);
	extractkVelSingleProbe(k3Earth, k3jw, NUMPERTURBED, k3Vel);

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
	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
//	std::cout << std::endl << "k4" << std::endl;
	singleProbeRK4Step(&k4jw, &k4Earth, NUMPERTURBED, k3Accel, k3Vel, dt);
//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;
	extractkAccelSingleProbe(k4Earth, sun, k4jw, NUMPERTURBED, k4Accel);

//	std::cout << std::endl << std::endl << "THE ULTIMATE DEBUG" << std::endl;
//	std::cout << std::endl << "DEBUGGING K4 AFTER STEP" << std::endl;
//	std::cout << "Printing k4 object positions" << std::endl;
//	std::cout << "k4Earth: x: " << k4Earth.pos.x << " " << k4Earth.pos.y << std::endl;
//	std::cout << "k4Earth: v: " << k4Earth.vel.x << " " << k4Earth.vel.y << std::endl;
//	std::cout << "k4JamesWebb: x: " << k4jw.pos.x << " " << k4jw.pos.y << std::endl;
//	std::cout << "k4JamesWebb: v: " << k4jw.vel.x << " " << k4jw.vel.y << std::endl << std::endl << std::endl << std::endl;
	extractkVelSingleProbe(k4Earth, k4jw, NUMPERTURBED, k4Vel);

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
	// Fill up weighted k's
	for (int i = 0; i < NUMPERTURBED + 4; ++i)
	{
		weightedKAccel[i] = {k1Accel[i].x + 2 * k2Accel[i].x + 2 * k3Accel[i].x + k4Accel[i].x, k1Accel[i].y + 2 * k2Accel[i].y + 2 * k3Accel[i].y + k4Accel[i].y};
		weightedKVel[i] = {k1Vel[i].x + 2 * k2Vel[i].x + 2 * k3Vel[i].x + k4Vel[i].x, k1Vel[i].y + 2 * k2Vel[i].y + 2 * k3Vel[i].y + k4Vel[i].y};
	}
	// Step real objects forward dt / 6 using weighted k's
	singleProbeRK4Step(jw, earth, NUMPERTURBED, weightedKAccel, weightedKVel, dt / 6.0);

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
	free(k2jw.perturbed); free(k2jw.h2); free(k3jw.perturbed); free(k3jw.h2); free(k4jw.perturbed); free(k4jw.h2);
	free(k1Accel); free(k2Accel); free(k3Accel); free(k4Accel); free(weightedKAccel); free(k1Vel); free(k2Vel); free(k3Vel); free(k4Vel); free(weightedKVel);

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
	for (int i = 0; i < NUMPERTURBED; ++i) { copyProbe(&jw->perturbed[i], *jw); }

	// Read from file
	// Desktop path
	FILE *f = fopen("D:\\lorenz3dChaos\\unitSphere.csv", "r");
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
