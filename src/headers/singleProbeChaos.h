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
void singleProbeRK4(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, double dt);
// Calculates and stores gravitational accelerations on earth, james webb, perturbations, and h2
void calculateAccels(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels);
// Stores velocities in kVels
void extractVel(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthVel, vector2 *jwVel, vector2 *perturbedVels, vector2 *h2Vels);
// Calculates and stores gravitational accelerations on earth, james webb, perturbations, and h2. Helper function for RK4
void extractkAccelSingleProbe(Planet earth, Planet sun, Probe jw, int NUMPERTURBED, vector2 *accels);
// Stores velocities in kVel: Helper function for RK4
void extractkVelSingleProbe(Planet earth, Probe jw, int NUMPERTURBED, vector2 *vels);
// Moves Objects forward using kAccel and velocity
void rk4Step(Planet *earth, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels, double dt);
// Moves objects forward using kAccel and kVel: Helper function for RK4
void singleProbeRK4Step(Probe *jwCopy, Planet *earthCopy, int NUMPERTURBED, vector2 *accels, vector2 *vels, double dt);
// Handles the timestep processes for calculating h1 and h2 lyapunov exponents
void lyapunovChaosStuff(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea, double RADIUS);
// Gram schmidt to re-orthonormalize unit vectors
void gramSchmidtRenormalizationh2(Probe *jw, int NUMPERTURBED, vector2 *v1, vector2 *v2, double RADIUS);
// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovs(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t);



#endif /* HEADERS_SINGLEPROBECHAOS_H_ */
