/*
 * unitTests.h
 *
 *  Created on: Feb 19, 2023
 *      Author: norms
 */

#ifndef HEADERS_UNITTESTS_H_
#define HEADERS_UNITTESTS_H_

#include "vector2.h"
#include "singleProbeChaos.h"
#include "PlanetAndProbe.h"

#include <iomanip>

// Multi Probe Chaos Unit Test Functions

// Prints out the perturbed offsets from the unit circle generation
void perturbedOffsetTest(vector2 *perturbedOffsets, int NUMPERTURBED);

// Prints out the James Webb positions that populate the grid with specified width and height
void probePositionsTest(vector2 *probePositions, int numProbes, int width, int height);

// Unit test for function that initializes probes, their perturbed probes and h2 probes.
void initializeProbesTest(Probe *probes, vector2 *probePositions, vector2 *probeVels, vector2 *perturbedOffsets, int numProbes, int NUMPERTURBED, double RADIUS);

// Returns a vector2 array of probe velocities such that they all have the same radial velocity
vector2 *velocitySpectrumTest(vector2 *probePositions, Planet *sun, int numProbes);

// Stores velocities in kVels
void multiProbeExtractVelsTest(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthVel, vector2 *probeVels, vector2 **perturbedVels, vector2 **h2Vels);
// Calculates and stores gravitational accelerations on earth, multiple probes, perturbations, and h2: RK4 Function
void multiProbeCalcAccelsTest(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels);
// Moves Objects forward using kAccel and velocity
void multiProbeRK4StepTest(Planet *earth, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels, double dt);

// Handles lyapunov calculations at every timestep to calculate h1 and h2
void multiProbeLyapunovStuffTest(Probe *probes, double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double RADIUS, double initialArea);

// Extracts the lyapunov results from a multi probe simulation
void multiProbeExtractLyapunovsTest(double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double **lyapunovs, double t);



// Single Probe Chaos
// Lyapunov unit tests
void gramSchmidtRenormalizationh2Test(Probe *jw, int NUMPERTURBED, vector2 *v1, vector2 *v2, double RADIUS);
void lyapunovChaosStuffTest(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea, double RADIUS);
// RK4 unit tests
void rk4StepTest(Planet *earth, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels, double dt);
void extractVelTest(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthVel, vector2 *jwVel, vector2 *perturbedVels, vector2 *h2Vels);
void calcAccelTest(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels);
// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovsTest(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t);


#endif /* HEADERS_UNITTESTS_H_ */
