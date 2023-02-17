/*
 * singleProbeChaos.h
 *
 *  Created on: Feb 14, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEPROBECHAOS_H_
#define HEADERS_SINGLEPROBECHAOS_H_

#include "PlanetAndProbe.h"

// Main function, extract two lyapunov exponents given a single probe with perturbations
double *singleProbeChaos(Planet *earth, Planet sun, Probe *jw, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, double RADIUS);

// Functions utilized by main function

// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid. Reads from file
void initializePerturbations(Probe *jw, int NUMPERTURBED);

// RK4 for sun earth system with single james webb with perturbations
void singleProbeRK4(Planet *earth, Planet sun, Probe *jw, int NUMPERTURBED, double dt);
// Calculates and stores gravitational accelerations on earth, james webb, perturbations, and h2. Helper function for RK4
void extractkAccelSingleProbe(Planet earth, Planet sun, Probe jw, int NUMPERTURBED, vector2 *accels);
// Stores velocities in kVel: Helper function for RK4
void extractkVelSingleProbe(Planet earth, Probe jw, int NUMPERTURBED, vector2 *vels);
// Moves objects forward using kAccel and kVel: Helper function for RK4
void singleProbeRK4Step(Probe *jwCopy, Planet *earthCopy, int NUMPERTURBED, vector2 *accels, vector2 *vels, double dt);

// Handles the timestep processes for calculating h1 and h2 lyapunov exponents
void lyapunovChaosStuffSingle(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea);
// Gram schmidt to re-orthonormalize unit vectors
void gramSchmidtRenormalizationh2Single(Probe *jw, int NUMPERTURBED, vector2 v1, vector2 v2);
// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovSingleProbe(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t);
// Writes position data to csv file for all objects in scene for plotting
void dataToFile(vector2 **data, int dataLen, int NUMPERTURBED);



#endif /* HEADERS_SINGLEPROBECHAOS_H_ */
