/*
 * multiProbeChaos.cpp
 *
 *  Created on: Feb 20, 2023
 *      Author: norms
 */


#include "../headers/multiProbeChaos.h"

// Main function, extracts lyapunov exponents for probes spanning a grid of specified resolution
double **multiProbeChaos(Planet *earth, Planet *sun, vector2 topLeft, vector2 bottomRight, int width, int height, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, double RADIUS)
{
	// Number of probes
	int numProbes = width * height;

	// Allocate space for lyapunov exponent 2d array
	double **lyapunovs = (double **)malloc(numProbes * sizeof(double *));

	// h1 lyapunov calculation: Orbital separation
	double **LfSums = (double **)malloc(numProbes * sizeof(double *));

	// h2 lyapunov Calculation: Sum of area
	double *areaSums = (double *)malloc(numProbes * sizeof(double));

	for (int i = 0; i < numProbes; ++i)
	{
		// Allocate space for two lyapunov exponents per probe
		lyapunovs[i] = (double *)malloc(2 * sizeof(double));

		// Allocate space for LfSums
		LfSums[i] = (double *)malloc(NUMPERTURBED * sizeof(double));

		// Initialize areasum to zero
		areaSums[i] = 0;

		// Initialize all LfSums to zero
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			LfSums[i][j] = 0;
		}
	}

	// Write the circle of points to file
	writeCircleToFile(NUMPERTURBED, RADIUS);

	// Read the circle into a data array
	vector2 *perturbedOffsets = returnPerturbationsCircle(NUMPERTURBED);

	// Create array of probe positions populating the grid space
	vector2 *probePositions = partitionPositions(topLeft, bottomRight, width, height);

	// Prints probe positions
	probePositionsTest(probePositions, numProbes, width, height);

	// Create array of probe velocities
	vector2 *probeVels = velocitySpectrum(probePositions, sun, numProbes);

	// Allocate memory for all James Webb Telescope Probes
	Probe *probes = (Probe *)malloc(numProbes * sizeof(Probe));

	// Initialize the James Webb Probes, their perturbed probes, and h2 probes
	initializeProbes(probes, probePositions, probeVels, perturbedOffsets, numProbes, NUMPERTURBED, RADIUS);

	// Store data: Store all data of probe and perturbations every ITERATIONS step: Useful for plotting.
	// Number of data buckets needed
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;

	// Position Data for earth
	vector2 *earthData = (vector2 *)malloc(n * sizeof(vector2));

	// Position Data for james webb probes
	vector2 **probesData = (vector2 **)malloc(n * sizeof(vector2 *));

	// Data for perturbations for all probes: 3d array
	vector2 ***perturbedData = (vector2 ***)malloc(n * sizeof(vector2 **));

	// Data for h2 for all probes: 3d array
	vector2 ***h2Data = (vector2 ***)malloc(n * sizeof(vector2 **));

	// Allocate memory for all data arrays
	for (int i = 0; i < n; ++i)
	{
		probesData[i] = (vector2 *)malloc(numProbes * sizeof(vector2));
		perturbedData[i] = (vector2 **)malloc(numProbes * sizeof(vector2 *));
		h2Data[i] = (vector2 **)malloc(numProbes * sizeof(vector2 *));

		for (int j = 0; j < numProbes; ++j)
		{
			perturbedData[i][j] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
			h2Data[i][j] = (vector2 *)malloc(2 * sizeof(vector2));
		}
	}

	// Loop
	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{
		// Store data
		// Store Data in different arrays
		if (c % ITERATIONS == 0)
		{
			std::cout << "c: " << c << std::endl;
			// Earth
			earthData[dc] = {earth->pos.x, earth->pos.y};

			// James Webb Probes:
			for (int i = 0; i < numProbes; ++i)
			{
				// Position
				probesData[dc][i] = {probes[i].pos.x, probes[i].pos.y};

				// Perturbed probes:
				for (int j = 0; j < NUMPERTURBED; ++j)
				{
					// Position
					perturbedData[dc][i][j] = {probes[i].perturbed[j].pos.x, probes[i].perturbed[j].pos.y};
				}

				// h2 probes:
				for (int j = 0; j < 2; ++j)
				{
					// Position
					h2Data[dc][i][j] = {probes[i].h2[j].pos.x, probes[i].h2[j].pos.y};
				}
			}

			// Increment the data count
			++dc;
		}

		// Numerical Integrator: RK4
		multiProbeRK4(earth, sun, probes, numProbes, NUMPERTURBED, dt);

		// Lyapunov exponent timestep calculations
//		multiProbeLyapunovStuff(probes, LfSums, areaSums, numProbes, NUMPERTURBED, RADIUS, RADIUS * RADIUS);

		// Increment c
		++c;
	}

	// Collect the last bucket of data
	// Earth
	earthData[n - 1] = {earth->pos.x, earth->pos.y};

	// James Webb Probes:
	for (int i = 0; i < numProbes; ++i)
	{
		// Position
		probesData[n - 1][i] = {probes[i].pos.x, probes[i].pos.y};

		// Perturbed probes:
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Position
			perturbedData[n - 1][i][j] = {probes[i].perturbed[j].pos.x, probes[i].perturbed[j].pos.y};
		}

		// h2 probes:
		for (int j = 0; j < 2; ++j)
		{
			// Position
			h2Data[n - 1][i][j] = {probes[i].h2[j].pos.x, probes[i].h2[j].pos.y};
		}
	}


	// Extract the lyapunovs
	// multiProbeExtractLyapunovs(LfSums, areaSums, numProbes, NUMPERTURBED, dt * (double)totalSteps);
	multiProbeExtractLyapunovs(LfSums, areaSums, numProbes, NUMPERTURBED, lyapunovs, dt * (double)totalSteps);

	// Write data to file
	// earthAndProbesToFile(earthData, probesData, n, numProbes);
	// lyapunovsToFile(lyapunovs, width, height);

	multiProbePerturbationEvolution(probesData, perturbedData, n, numProbes, NUMPERTURBED);

	// earthAndProbesToFile(earthData, probesData, perturbedData, h2Data, n, numProbes, NUMPERTURBED);


	// Return the lyapunovs
	return lyapunovs;

}

// Extracts the lyapunov results from a multi probe simulation
void multiProbeExtractLyapunovs(double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double **lyapunovs, double t)
{
	// Loop through probes
	for (int i = 0; i < numProbes; ++i)
	{
		// Calculate and store h1: max lyapunov

		// Find maximum sum of all perturbed probes
		double max = LfSums[i][0];

		for (int j = 1; j < NUMPERTURBED; ++j)
		{
			if (LfSums[i][j] > max)
			{
				max = LfSums[i][j];
			}
		}

		// The maximum lyapunov exponent is the maximum divided by the total elapsed time
		lyapunovs[i][0] = max / t;

		// Calculate and store the second lyapunov exponent
		lyapunovs[i][1] = (areaSums[i] / t) - lyapunovs[i][0];
	}
}
// Handles lyapunov calculations at every timestep to calculate h1 and h2
void multiProbeLyapunovStuff(Probe *probes, double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double RADIUS, double initialArea)
{
	// h1: Max lyapunov

	// Loop through all probes and perturbed Probes
	for (int i = 0; i < numProbes; ++i)
	{
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Calculate Lf of each perturbation, and renormalize the position vector
			// Vector pointing from probe to perturbed probe
			vector2 v = {probes[i].perturbed[j].pos.x - probes[i].pos.x, probes[i].perturbed[j].pos.y - probes[i].pos.y};

			// Find the distance of the vector and scale it relative to the inital distance
			double dist = mag(v) / RADIUS;

			// Take the natural log
			double logDist = log(dist);

			// Add to the running sum for this perturbed probe
			LfSums[i][j] += logDist;

			// Renormalize the perturbed probe preserving the direction vector
			vector2 uv = unitVector(v);

			// Scale unit vector by the initial distance
			uv.x *= RADIUS;
			uv.y *= RADIUS;

			// Renormalize the perturbed probe
			probes[i].perturbed[j].pos.x = probes[i].pos.x + uv.x; probes[i].perturbed[j].pos.y = probes[i].pos.y + uv.y;
		}

		// h2 lyapunov calculation:

		// Two vectors pointing from probe to each h2 probe
		vector2 v1 = {probes[i].h2[0].pos.x - probes[i].pos.x, probes[i].h2[0].pos.y - probes[i].pos.y};
		vector2 v2 = {probes[i].h2[1].pos.x - probes[i].pos.x, probes[i].h2[1].pos.y - probes[i].pos.y};

		// Calculate area of parallelogram formed by two vectors: A = ||a|| ||b|| sin(ab)
		double area = mag(v1) * mag(v2) * sin(angle(v1, v2));

		// Scale the area relative to the initial area
		double scaledArea = area / initialArea;

		// Take natural log
		double logArea = log(scaledArea);

		// Add to the areaSum
		areaSums[i] += logArea;

		// Renormalize the vectors using Gram Schmidt Algorithm to preserve direction:
		// First orthonormal basis vector
		vector2 e1 = unitVector(v1);

		// Dot product between v2 and e1
		double v2Dote1 = dotProd(v2, e1);

		// Calculate the orthogonal basis vector
		vector2 y2 = {v2.x - (v2Dote1 * e1.x), v2.y - (v2Dote1 * e1.y)};

		// Take unit vector to get second orthonormal basis vector
		vector2 e2 = unitVector(y2);

		// Scale the basis vectors by the initial distance
		e1.x *= RADIUS;
		e1.y *= RADIUS;
		e2.x *= RADIUS;
		e2.y *= RADIUS;

		// Renormalize the h2 probes
		probes[i].h2[0].pos = {probes[i].pos.x + e1.x, probes[i].pos.y + e1.y};
		probes[i].h2[1].pos = {probes[i].pos.x + e2.x, probes[i].pos.y + e2.y};
	}
}

// RK4 for sun earth system with multiple James Webb Telescopes with perturbed probes, and h2 probes
void multiProbeRK4(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, double dt)
{
	// Initialize earth Copies: k2 - k4
	Planet k2Earth, k3Earth, k4Earth;

	// Allocate memory for James Webb Array Copies: k2 - k4
	Probe *k2Probes = (Probe *)malloc(numProbes * sizeof(Probe));
	Probe *k3Probes = (Probe *)malloc(numProbes * sizeof(Probe));
	Probe *k4Probes = (Probe *)malloc(numProbes * sizeof(Probe));

	// Initialize kAccels and kVels:
	// Earth:
	vector2 k1EarthAccel, k2EarthAccel, k3EarthAccel, k4EarthAccel, weightedEarthAccel;
	vector2 k1EarthVel, k2EarthVel, k3EarthVel, k4EarthVel, weightedEarthVel;

	// James Webb Probes:
	vector2 *k1ProbeAccels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k2ProbeAccels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k3ProbeAccels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k4ProbeAccels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *weightedProbeAccels = (vector2 *)malloc(numProbes * sizeof(vector2));

	vector2 *k1ProbeVels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k2ProbeVels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k3ProbeVels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *k4ProbeVels = (vector2 *)malloc(numProbes * sizeof(vector2));
	vector2 *weightedProbeVels = (vector2 *)malloc(numProbes * sizeof(vector2));

	// Perturbations: Each are a 2d array
	vector2 **k1PerturbedAccels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k2PerturbedAccels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k3PerturbedAccels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k4PerturbedAccels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **weightedPerturbedAccels = (vector2 **)malloc(numProbes * sizeof(vector2 *));

	vector2 **k1PerturbedVels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k2PerturbedVels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k3PerturbedVels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k4PerturbedVels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **weightedPerturbedVels = (vector2 **)malloc(numProbes * sizeof(vector2 *));

	// h2 probes: Each are a 2d array
	vector2 **k1h2Accels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k2h2Accels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k3h2Accels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k4h2Accels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **weightedh2Accels = (vector2 **)malloc(numProbes * sizeof(vector2 *));

	vector2 **k1h2Vels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k2h2Vels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k3h2Vels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **k4h2Vels = (vector2 **)malloc(numProbes * sizeof(vector2 *));
	vector2 **weightedh2Vels = (vector2 **)malloc(numProbes * sizeof(vector2 *));

	// Copy data into new copies
	// Earth:
	copyPlanet(&k2Earth, *earth);
	copyPlanet(&k3Earth, *earth);
	copyPlanet(&k4Earth, *earth);

	// Probes:
	for (int i = 0; i < numProbes; ++i)
	{
		// Copy data into James Webb copies
		copyProbe(&k2Probes[i], probes[i]);
		copyProbe(&k3Probes[i], probes[i]);
		copyProbe(&k4Probes[i], probes[i]);

		// Allocate memory for perturbed probes
		k2Probes[i].perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
		k3Probes[i].perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));
		k4Probes[i].perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));

		// Copy data into perturbed probes
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			copyProbe(&k2Probes[i].perturbed[j], probes[i].perturbed[j]);
			copyProbe(&k3Probes[i].perturbed[j], probes[i].perturbed[j]);
			copyProbe(&k4Probes[i].perturbed[j], probes[i].perturbed[j]);
		}

		// Allocate memory for h2 probes
		k2Probes[i].h2 = (Probe *)malloc(2 * sizeof(Probe));
		k3Probes[i].h2 = (Probe *)malloc(2 * sizeof(Probe));
		k4Probes[i].h2 = (Probe *)malloc(2 * sizeof(Probe));

		// Copy data into h2 probes
		for (int j = 0; j < 2; ++j)
		{
			copyProbe(&k2Probes[i].h2[j], probes[i].h2[j]);
			copyProbe(&k3Probes[i].h2[j], probes[i].h2[j]);
			copyProbe(&k4Probes[i].h2[j], probes[i].h2[j]);
		}

		// Allocate memory for:
		// Perturbed Accelerations
		k1PerturbedAccels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k2PerturbedAccels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k3PerturbedAccels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k4PerturbedAccels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		weightedPerturbedAccels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));

		// Perturbed Velocities
		k1PerturbedVels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k2PerturbedVels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k3PerturbedVels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		k4PerturbedVels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));
		weightedPerturbedVels[i] = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));

		// h2 Accelerations
		k1h2Accels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k2h2Accels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k3h2Accels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k4h2Accels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		weightedh2Accels[i] = (vector2 *)malloc(2 * sizeof(vector2));

		// h2 Velocities
		k1h2Vels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k2h2Vels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k3h2Vels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		k4h2Vels[i] = (vector2 *)malloc(2 * sizeof(vector2));
		weightedh2Vels[i] = (vector2 *)malloc(2 * sizeof(vector2));
	}

	// k1 Step:

//	std::cout << "k1 Step: " << std::endl << std::endl;

	// Calculate accelerations for all objects and store in k1Accels
	multiProbeCalcAccels(earth, sun, probes, numProbes, NUMPERTURBED, &k1EarthAccel, k1ProbeAccels, k1PerturbedAccels, k1h2Accels);

	// Extract velocities and store in k1Vels
	multiProbeExtractVels(earth, sun, probes, numProbes, NUMPERTURBED, &k1EarthVel, k1ProbeVels, k1PerturbedVels, k1h2Vels);

	// k2 Step:

//	std::cout << "k2 Step: " << std::endl << std::endl;

	// Step k2 objects forward by dt / 2 using k1
	multiProbeRK4Step(&k2Earth, k2Probes, numProbes, NUMPERTURBED, &k1EarthAccel, k1ProbeAccels, k1PerturbedAccels, k1h2Accels, dt / 2.0);

	// Calculate accelerations for all objects and store in k2Accels
	multiProbeCalcAccels(&k2Earth, sun, k2Probes, numProbes, NUMPERTURBED, &k2EarthAccel, k2ProbeAccels, k2PerturbedAccels, k2h2Accels);

	// Extract velocities and store in k2Vels
	multiProbeExtractVels(&k2Earth, sun, k2Probes, numProbes, NUMPERTURBED, &k2EarthVel, k2ProbeVels, k2PerturbedVels, k2h2Vels);

	// k3 Step:

//	std::cout << "k3 Step: " << std::endl << std::endl;

	// Step k3 objects forward by dt / 2 using k2
	multiProbeRK4Step(&k3Earth, k3Probes, numProbes, NUMPERTURBED, &k2EarthAccel, k2ProbeAccels, k2PerturbedAccels, k2h2Accels, dt / 2.0);

	// Calculate accelerations for all objects and store in k3Accels
	multiProbeCalcAccels(&k3Earth, sun, k3Probes, numProbes, NUMPERTURBED, &k3EarthAccel, k3ProbeAccels, k3PerturbedAccels, k3h2Accels);

	// Extract velocities and store in k3Vels
	multiProbeExtractVels(&k3Earth, sun, k3Probes, numProbes, NUMPERTURBED, &k3EarthVel, k3ProbeVels, k3PerturbedVels, k3h2Vels);

	// k4 Step:

//	std::cout << "k4 Step: " << std::endl << std::endl;

	// Step k4 objects forward by dt using k3
	multiProbeRK4Step(&k4Earth, k4Probes, numProbes, NUMPERTURBED, &k3EarthAccel, k3ProbeAccels, k3PerturbedAccels, k3h2Accels, dt);

	// Calculate accelerations for all objects and store in k4Accels
	multiProbeCalcAccels(&k4Earth, sun, k4Probes, numProbes, NUMPERTURBED, &k4EarthAccel, k4ProbeAccels, k4PerturbedAccels, k4h2Accels);

	// Extract velocities and store in k4Vels
	multiProbeExtractVels(&k4Earth, sun, k4Probes, numProbes, NUMPERTURBED, &k4EarthVel, k4ProbeVels, k4PerturbedVels, k4h2Vels);

	// RK4 Step:
//	std::cout << "RK4 Step" << std::endl << std::endl;
	// Calculate and store weighted accelerations and velocities

	// Earth:
	// Acceleration
	weightedEarthAccel.x = k1EarthAccel.x + 2 * k2EarthAccel.x + 2 * k3EarthAccel.x + k4EarthAccel.x;
	weightedEarthAccel.y = k1EarthAccel.y + 2 * k2EarthAccel.y + 2 * k3EarthAccel.y + k4EarthAccel.y;

//	std::cout << "Earth" << std::endl;
//	std::cout << "Weighted acceleration: x: " << weightedEarthAccel.x << " y: " << weightedEarthAccel.y << std::endl;

	// Velocity
	weightedEarthVel.x = k1EarthVel.x + 2 * k2EarthVel.x + 2 * k3EarthVel.x + k4EarthVel.x;
	weightedEarthVel.y = k1EarthVel.y + 2 * k2EarthVel.y + 2 * k3EarthVel.y + k4EarthVel.y;

//	std::cout << "Weighted velocity: x: " << weightedEarthVel.x << " y: " << weightedEarthVel.y << std::endl;

	// Probes:

//	std::cout << "Probes:" << std::endl << std::endl;

	for (int i = 0; i < numProbes; ++i)
	{
		// Acceleration
		weightedProbeAccels[i].x = k1ProbeAccels[i].x + 2 * k2ProbeAccels[i].x + 2 * k3ProbeAccels[i].x + k4ProbeAccels[i].x;
		weightedProbeAccels[i].y = k1ProbeAccels[i].y + 2 * k2ProbeAccels[i].y + 2 * k3ProbeAccels[i].y + k4ProbeAccels[i].y;

//		std::cout << "Probe: " << i << std::endl;
//		std::cout << "Weighted acceleration: x: " << weightedProbeAccels[i].x << " y: " << weightedProbeAccels[i].y << std::endl;
		// Velocity
		weightedProbeVels[i].x = k1ProbeVels[i].x + 2 * k2ProbeVels[i].x + 2 * k3ProbeVels[i].x + k4ProbeVels[i].x;
		weightedProbeVels[i].y = k1ProbeVels[i].y + 2 * k2ProbeVels[i].y + 2 * k3ProbeVels[i].y + k4ProbeVels[i].y;

		// Perturbed Probes:
//		std::cout << "Perturbed Probes:" << std::endl << std::endl;
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Acceleration
			weightedPerturbedAccels[i][j].x = k1PerturbedAccels[i][j].x + 2 * k2PerturbedAccels[i][j].x + 2 * k3PerturbedAccels[i][j].x + k4PerturbedAccels[i][j].x;
			weightedPerturbedAccels[i][j].y = k1PerturbedAccels[i][j].y + 2 * k2PerturbedAccels[i][j].y + 2 * k3PerturbedAccels[i][j].y + k4PerturbedAccels[i][j].y;

//			std::cout << "Perturbed Probe: " << j << std::endl;
//			std::cout << "Weighted acceleration: x: " << weightedPerturbedAccels[i][j].x << " y: " << weightedPerturbedAccels[i][j].y << std::endl;

			// Velocity
			weightedPerturbedVels[i][j].x = k1PerturbedVels[i][j].x + 2 * k2PerturbedVels[i][j].x + 2 * k3PerturbedVels[i][j].x + k4PerturbedVels[i][j].x;
			weightedPerturbedVels[i][j].y = k1PerturbedVels[i][j].y + 2 * k2PerturbedVels[i][j].y + 2 * k3PerturbedVels[i][j].y + k4PerturbedVels[i][j].y;

//			std::cout << "Weighted velocity: x: " << weightedPerturbedVels[i][j].x << " y: " << weightedPerturbedVels[i][j].y << std::endl << std::endl;
		}

		// h2 Probes:
//		std::cout << "h2 Probes:" << std::endl << std::endl;
		for (int j = 0; j < 2; ++j)
		{
			// Acceleration
			weightedh2Accels[i][j].x = k1h2Accels[i][j].x + 2 * k2h2Accels[i][j].x + 2 * k3h2Accels[i][j].x + k4h2Accels[i][j].x;
			weightedh2Accels[i][j].y = k1h2Accels[i][j].y + 2 * k2h2Accels[i][j].y + 2 * k3h2Accels[i][j].y + k4h2Accels[i][j].y;

//			std::cout << "h2 Probe: " << j << std::endl;
//			std::cout << "Weighted acceleration: x: " << weightedh2Accels[i][j].x << " y: " << weightedh2Accels[i][j].y << std::endl;

			// Velocity
			weightedh2Vels[i][j].x = k1h2Vels[i][j].x + 2 * k2h2Vels[i][j].x + 2 * k3h2Vels[i][j].x + k4h2Vels[i][j].x;
			weightedh2Vels[i][j].y = k1h2Vels[i][j].y + 2 * k2h2Vels[i][j].y + 2 * k3h2Vels[i][j].y + k4h2Vels[i][j].y;

//			std::cout << "Weighted velocity: x: " << weightedh2Vels[i][j].x << " y: " << weightedh2Vels[i][j].y << std::endl;
		}
	}
//	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;


	// Step real objects forward by dt / 6.0
//	std::cout << "Stepping objects forward by weighted accels and velocities" << std::endl << std::endl;

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

	// Probes:
//	std::cout << "Probes:" << std::endl << std::endl;
	for (int i = 0; i < numProbes; ++i)
	{
		// Velocity

//		std::cout << "Probe: " << i << std::endl;
//		std::cout << "Old velocities: vx: " << probes[i].vel.x << " vy: " << probes[i].vel.y << std::endl;

		probes[i].vel.x += weightedProbeAccels[i].x * dt / 6.0;
		probes[i].vel.y += weightedProbeAccels[i].y * dt / 6.0;

		// Position

//		std::cout << "New velocities: vx: " << probes[i].vel.x << " vy: " << probes[i].vel.y << std::endl;
//		std::cout << "Old positions: x: " << probes[i].pos.x << " y: " << probes[i].pos.y << std::endl;

		probes[i].pos.x += weightedProbeVels[i].x * dt / 6.0;
		probes[i].pos.y += weightedProbeVels[i].y * dt / 6.0;

//		std::cout << "New positions: x: " << probes[i].pos.x << " y: " << probes[i].pos.y << std::endl;

		// Perturbed Probes:
//		std::cout << "Perturbed Probes:" << std::endl << std::endl;
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Velocity
//			std::cout << "Perturbed Probe: " << j << std::endl;
//			std::cout << "Old velocities: vx: " << probes[i].perturbed[j].vel.x << " vy: " << probes[i].perturbed[j].vel.y << std::endl;

			probes[i].perturbed[j].vel.x += weightedPerturbedAccels[i][j].x * dt / 6.0;
			probes[i].perturbed[j].vel.y += weightedPerturbedAccels[i][j].y * dt / 6.0;

			// Position
//			std::cout << "New velocities: vx: " << probes[i].perturbed[j].vel.x << " vy: " << probes[i].perturbed[j].vel.y << std::endl;
//			std::cout << "Old positions: x: " << probes[i].perturbed[j].pos.x << " y: " << probes[i].perturbed[j].pos.y << std::endl;

			probes[i].perturbed[j].pos.x += weightedPerturbedVels[i][j].x * dt / 6.0;
			probes[i].perturbed[j].pos.y += weightedPerturbedVels[i][j].y * dt / 6.0;

//			std::cout << "New positions: x: " << probes[i].perturbed[j].pos.x << " y: " << probes[i].perturbed[j].pos.y << std::endl;
		}

		// h2 Probes:
//		std::cout << "h2 Probes:" << std::endl << std::endl;
		for (int j = 0; j < 2; ++j)
		{
			// Velocity
//			std::cout << "h2 Probe: " << j << std::endl;
//			std::cout << "Old velocities: vx: " << probes[i].h2[j].vel.x << " vy: " << probes[i].h2[j].vel.y << std::endl;

			probes[i].h2[j].vel.x += weightedh2Accels[i][j].x * dt / 6.0;
			probes[i].h2[j].vel.y += weightedh2Accels[i][j].y * dt / 6.0;

			// Position
//			std::cout << "New velocities: vx: " << probes[i].h2[j].vel.x << " vy: " << probes[i].h2[j].vel.y << std::endl;
//			std::cout << "Old positions: x: " << probes[i].h2[j].pos.x << " y: " << probes[i].h2[j].pos.y << std::endl;

			probes[i].h2[j].pos.x += weightedh2Vels[i][j].x * dt / 6.0;
			probes[i].h2[j].pos.y += weightedh2Vels[i][j].y * dt / 6.0;

//			std::cout << "New positions: x: " << probes[i].h2[j].pos.x << " y: " << probes[i].h2[j].pos.y << std::endl;
		}
//		std::cout << std::endl;
	}

	// Free all of the dynamically allocated arrays:

	// The 2d arrays

	for (int i = 0; i < numProbes; ++i)
	{
		// Free perturbed probes
		free(k2Probes[i].perturbed);
		free(k3Probes[i].perturbed);
		free(k4Probes[i].perturbed);

		// Free h2 probes
		free(k2Probes[i].h2);
		free(k3Probes[i].h2);
		free(k4Probes[i].h2);

		// Perturbed:

		// Accelerations
		free(k1PerturbedAccels[i]);
		free(k2PerturbedAccels[i]);
		free(k3PerturbedAccels[i]);
		free(k4PerturbedAccels[i]);
		free(weightedPerturbedAccels[i]);

		// Velocities
		free(k1PerturbedVels[i]);
		free(k2PerturbedVels[i]);
		free(k3PerturbedVels[i]);
		free(k4PerturbedVels[i]);
		free(weightedPerturbedVels[i]);

		// h2:

		// Accelerations
		free(k1h2Accels[i]);
		free(k2h2Accels[i]);
		free(k3h2Accels[i]);
		free(k4h2Accels[i]);
		free(weightedh2Accels[i]);

		// Velocities
		free(k1h2Vels[i]);
		free(k2h2Vels[i]);
		free(k3h2Vels[i]);
		free(k4h2Vels[i]);
		free(weightedh2Vels[i]);
	}

	// Probes:
	free(k2Probes);
	free(k3Probes);
	free(k4Probes);

	// Accelerations
	free(k1ProbeAccels);
	free(k2ProbeAccels);
	free(k3ProbeAccels);
	free(k4ProbeAccels);
	free(weightedProbeAccels);

	// Velocities
	free(k1ProbeVels);
	free(k2ProbeVels);
	free(k3ProbeVels);
	free(k4ProbeVels);
	free(weightedProbeVels);

	// Perturbed:

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

// Moves Objects forward using kAccel and velocity
void multiProbeRK4Step(Planet *earth, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels, double dt)
{
	// Push velocities forward using k Accel, and then push positions forward using current velocity

	// Earth:
	// Velocities
	earth->vel.x += earthAccel->x * dt;
	earth->vel.y += earthAccel->y * dt;

	// Positions
	earth->pos.x += earth->vel.x * dt;
	earth->pos.y += earth->vel.y * dt;

	// Probes:
	for (int i = 0; i < numProbes; ++i)
	{
		// Velocities
		probes[i].vel.x += probeAccels[i].x * dt;
		probes[i].vel.y += probeAccels[i].y * dt;

		// Positions
		probes[i].pos.x += probes[i].vel.x * dt;
		probes[i].pos.y += probes[i].vel.y * dt;

		// Perturbed Probes:
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Velocities
			probes[i].perturbed[j].vel.x += perturbedAccels[i][j].x * dt;
			probes[i].perturbed[j].vel.y += perturbedAccels[i][j].y * dt;

			// Positions
			probes[i].perturbed[j].pos.x += probes[i].perturbed[j].vel.x * dt;
			probes[i].perturbed[j].pos.y += probes[i].perturbed[j].vel.y * dt;
		}

		// h2 Probes:
		for (int j = 0; j < 2; ++j)
		{
			// Velocities
			probes[i].h2[j].vel.x += h2Accels[i][j].x * dt;
			probes[i].h2[j].vel.y += h2Accels[i][j].y * dt;

			// Positions
			probes[i].h2[j].pos.x += probes[i].h2[j].vel.x * dt;
			probes[i].h2[j].pos.y += probes[i].h2[j].vel.y * dt;
		}
	}
}

// Calculates and stores gravitational accelerations on earth, multiple probes, perturbations, and h2
void multiProbeCalcAccels(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels)
{
	// Gravitational Constant
	double G = 6.6743e-11;

	// Sun on earth:
	vector2 r = {earth->pos.x - sun->pos.x, earth->pos.y - sun->pos.y}; 	// Vector pointing from sun to earth
	double magr = mag(r);														// Magnitude of that vector
	vector2 rUnit = unitVector(r);												// Direction of the vector
	double g = (-1 * G * sun->m) / (magr * magr);								// Newtons law of gravity to get force magnitude
	// Multiply by unit vector for direction
	earthAccel->x = g * rUnit.x;
	earthAccel->y = g * rUnit.y;

	// James Webb Probes:
	for (int i = 0; i < numProbes; ++i)
	{
		// Add up accelerations from Sun, and from Earth

		// Sun on James Webb Probe
		r = {probes[i].pos.x - sun->pos.x, probes[i].pos.y - sun->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * sun->m) / (magr * magr);
		probeAccels[i].x = g * rUnit.x;
		probeAccels[i].y = g * rUnit.y;

		// Earth on James Webb Probe
		r = {probes[i].pos.x - earth->pos.x, probes[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		// Add to accelerations from the Sun
		probeAccels[i].x += g * rUnit.x;
		probeAccels[i].y += g * rUnit.y;

		// Perturbed Probes:
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Add up accelerations from Sun, and from Earth

			// Sun on perturbation
			r = {probes[i].perturbed[j].pos.x - sun->pos.x, probes[i].perturbed[j].pos.y - sun->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * sun->m) / (magr * magr);
			perturbedAccels[i][j].x = g * rUnit.x;
			perturbedAccels[i][j].y = g * rUnit.y;

			// Earth on perturbation
			r = {probes[i].perturbed[j].pos.x - earth->pos.x, probes[i].perturbed[j].pos.y - earth->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * earth->m) / (magr * magr);
			perturbedAccels[i][j].x += g * rUnit.x;
			perturbedAccels[i][j].y += g * rUnit.y;
		}

		// h2 Probes:
		for (int j = 0; j < 2; ++j)
		{
			// Add up accelerations from Sun, and from Earth

			// Sun on h2 probe
			r = {probes[i].h2[j].pos.x - sun->pos.x, probes[i].h2[j].pos.y - sun->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * sun->m) / (magr * magr);
			h2Accels[i][j].x = g * rUnit.x;
			h2Accels[i][j].y = g * rUnit.y;

			// Earth on h2 probe
			r = {probes[i].h2[j].pos.x - earth->pos.x, probes[i].h2[j].pos.y - earth->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * earth->m) / (magr * magr);
			h2Accels[i][j].x += g * rUnit.x;
			h2Accels[i][j].y += g * rUnit.y;
		}
	}
}

// Stores velocities in kVels
void multiProbeExtractVels(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthVel, vector2 *probeVels, vector2 **perturbedVels, vector2 **h2Vels)
{
	// Earth
	earthVel->x = earth->vel.x;
	earthVel->y = earth->vel.y;

	// Probes:
	for (int i = 0; i < numProbes; ++i)
	{
		probeVels[i].x = probes[i].vel.x;
		probeVels[i].y = probes[i].vel.y;

		// Perturbed Probes:
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			perturbedVels[i][j].x = probes[i].perturbed[j].vel.x;
			perturbedVels[i][j].y = probes[i].perturbed[j].vel.y;
		}

		// h2 Probes:
		for (int j = 0; j < 2; ++j)
		{
			h2Vels[i][j].x = probes[i].h2[j].vel.x;
			h2Vels[i][j].y = probes[i].h2[j].vel.y;
		}
	}
}

// Initializes all James Webb probes, their perturbed probes, and their h2 probes.
void initializeProbes(Probe *probes, vector2 *probePositions, vector2 *probeVels, vector2 *perturbedOffsets, int numProbes, int NUMPERTURBED, double RADIUS)
{
	for (int i = 0; i < numProbes; ++i)
	{
		// James Webb:
		// Position
		probes[i].pos = {probePositions[i].x, probePositions[i].y};

		// Velocity
		probes[i].vel = {probeVels[i].x, probeVels[i].y};

		// Perturbed Probes:
		// Allocate space for perturbed probes
		probes[i].perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));

		// Copy information over
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Position: Add the perturbed offset to the James Webb position
			probes[i].perturbed[j].pos = {probes[i].pos.x + perturbedOffsets[j].x, probes[i].pos.y + perturbedOffsets[j].y};

			// Velocity:
			probes[i].perturbed[j].vel = {probes[i].vel.x, probes[i].vel.y};
		}

		// h2 Orthonormal Probes:

		// Allocate space
		probes[i].h2 = (Probe *)malloc(2 * sizeof(Probe));

		// Copy information over
		for (int j = 0; j < 2; ++j)
		{
			copyProbe(&probes[i].h2[j], probes[i]);
		}

		// Add orthonormal offsets
		probes[i].h2[0].pos.x += RADIUS;
		probes[i].h2[1].pos.y += RADIUS;
	}
}

// Returns a vector2 array of probe velocities such that they all have the same radial velocity
vector2 *velocitySpectrum(vector2 *probePositions, Planet *sun, int numProbes)
{
	// Allocate space for probe velocities
	vector2 *probeVels = (vector2 *)malloc(numProbes * sizeof(vector2));

	// Loop through all probes
	for (int i = 0; i < numProbes; ++i)
	{
		// Vector from the probe to the sun
		vector2 r = {sun->pos.x - probePositions[i].x, sun->pos.y - probePositions[i].y};

		// Unit vector pointing from probe to sun
		vector2 rUnit = unitVector(r);

		// Magnitude
		double dist = mag(r);

		// Radial speed of the earth: radialSpeed = v / earthDist
		double earthSpeed = 29787.7;
		double earthToSunDist = 149597870700;
		double radialSpeed = earthSpeed / earthToSunDist;

		// Find speed of probe: v = dist * radialSpeed
		double probeSpeed = radialSpeed * dist;

		// Perpendicular vector to unit vector in the counter clockwise direction
		vector2 rUnitOrbital = {-1 * rUnit.y, rUnit.x};

		// Calculate and store velocity vector
		probeVels[i] = {probeSpeed * rUnitOrbital.x, probeSpeed * rUnitOrbital.y};
	}

	// Return the probe velocity array
	return probeVels;
}
// Returns an array of x, y probe positions distributed across a square for all probes
vector2 *partitionPositions(vector2 topLeft, vector2 bottomRight, int width, int height)
{
	// Find the increment in x and y
	double xIncr = (bottomRight.x - topLeft.x) / (double)(width - 1);
	double yIncr = (topLeft.y - bottomRight.y) / (double)(height - 1);

	// Allocate memory for the probe positions
	vector2 *positions = (vector2 *)malloc(width * height * sizeof(vector2));

	for (int i = 0; i < width * height; ++i)
	{
		positions[i].x = topLeft.x + xIncr * (i % width);
		positions[i].y = bottomRight.y + yIncr * (i / width);
	}

	return positions;
}

// Returns an array structure of circular perturbations in xz plane read from file. Faster to read from array for use over multiple probes. Struct: data[NUMPOINTS][x, z]
vector2 *returnPerturbationsCircle(int NUMPERTURBED)
{
	// Allocate memory for the array
	vector2 *perturbed = (vector2 *)malloc(NUMPERTURBED * sizeof(vector2));

	// Read from file
	FILE *f = fopen("unitCircle.csv", "r");
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);

	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);

	// Read from buffer into struct
	int pointCount = 0;
	char arr[30];
	int arrCount = 0;
	int i = 0;

	while (pointCount < NUMPERTURBED)
	{
		if (buffer[i] == '\n')
		{
			arr[arrCount] = '\0';
			perturbed[pointCount].y = std::atof(arr);
			pointCount++;
			arrCount = 0;
		}
		else if (buffer[i] == ',')
		{
			arr[arrCount] = '\0';
			perturbed[pointCount].x = std::atof(arr);
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

	return perturbed;
}

