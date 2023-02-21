//============================================================================
// Name        : 2dSunEarthChaos.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <math.h>
#include <iomanip>


#include "headers/singleProbeChaos.h"
#include "headers/multiProbeChaos.h"

int main() {

	// Timestep: JWST crashes into Earth at 28506500 seconds
	double dt = 50.0;
	int totalSteps = 450000;

	// Data storage
	int ITERATIONS = 225;

	// Initialize Earth and Sun
	Planet sun; sun.pos = {149597870700, 0}; sun.vel = {0, 0}; sun.m = 1.989e30;
	Planet earth; earth.pos = {0, 0}; earth.vel = {0, 29787.7}; earth.m = 5.972e24;

	// For multi probe chaos

	// Specify the top left and bottom right points of grid
	vector2 topLeft = {-1.7e9, 5e8};
	vector2 bottomRight = {-1.3e9, -5e8};

	// Specify the grid resolution
	int width = 3;
	int height = 3;

	// Number of perturbed probes
	int NUMPERTURBED = 36;

	// Initial distance of perturbed probes
	double RADIUS = 1000.0;


//	// Call the function
//	double **multiLyapunovs = multiProbeChaos(&earth, &sun, topLeft, bottomRight, width, height, dt, totalSteps, ITERATIONS, NUMPERTURBED, RADIUS);
//
//	// Free lyapunovs
//	for (int i = 0; i < width * height; ++i)
//	{
//		free(multiLyapunovs[i]);
//	}
//
//	free(multiLyapunovs);

	// For single probe chaos
	Probe jw;
	// Reinitialize earth
	earth.pos = {0, 0}; earth.vel = {0, 29787.7}; earth.m = 5.972e24;
	jw.pos = {-1.5e9, 0};
	jw.vel = {0, 30086.377713733};

	// Main function, extract two lyapunov exponents given a single probe with perturbations
	double *lyapunovs = singleProbeChaos(&earth, sun, &jw, dt, totalSteps, ITERATIONS, NUMPERTURBED, RADIUS);
	std::cout << "h1: " << std::setprecision(10) << lyapunovs[0] << " h2: " << std::setprecision(10) << lyapunovs[1];
	free(lyapunovs);
	return 0;
}
