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

#include "headers/circleGenerator.h"
#include "headers/PlanetAndProbe.h"
#include "headers/vector2.h"
#include "headers/singleProbeChaos.h"


int main() {

	// Timestep: JWST crashes into Earth at 28506500 seconds
	double dt = 1.0;
	int totalSteps = 50;

	// Data storage
	int ITERATIONS = 1;

	// Initialize Planets and Probes
	Planet sun; sun.pos = {149597870700, 0}; sun.vel = {0, 0}; sun.m = 1.989e30;
	Planet earth; earth.pos = {0, 0}; earth.vel = {0, 29787.7}; earth.m = 5.972e24;
	Probe jw; jw.pos = {-1.5e9, 0}; jw.vel = {0, 30173.943};

	// Single Probe Main Function
	int NUMPERTURBED = 5; 					// Number of perturbations
	double RADIUS = 1.0;					// Distance of perturbations
	// Main function, extract two lyapunov exponents given a single probe with perturbations
	double *lyapunovs = singleProbeChaos(&earth, sun, &jw, dt, totalSteps, ITERATIONS, NUMPERTURBED, RADIUS);
	std::cout << "h1: " << std::setprecision(10) << lyapunovs[0] << " h2: " << std::setprecision(10) << lyapunovs[1];
	free(lyapunovs);
	std::cout << "Done" << std::endl;
	return 0;
}
