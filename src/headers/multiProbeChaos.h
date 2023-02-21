/*
 * multiProbeChaos.h
 *
 *  Created on: Feb 20, 2023
 *      Author: norms
 */

#ifndef HEADERS_MULTIPROBECHAOS_H_
#define HEADERS_MULTIPROBECHAOS_H_

#include "circleGenerator.h"
#include "vector2.h"
#include "PlanetAndProbe.h"
#include "toFile.h"
#include "unitTests.h"
#include "multiProbeToFile.h"

#include <iostream>
#include <iomanip>
#include "math.h"


// Main function, extracts lyapunov exponents for probes spanning a grid of specified resolution
double **multiProbeChaos(Planet *earth, Planet *sun, vector2 topLeft, vector2 bottomRight, int width, int height, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, double RADIUS);

// Returns an array of x, y probe positions distributed across a square for all probes
vector2 *partitionPositions(vector2 topLeft, vector2 bottomRight, int width, int height);
// Returns a vector2 array of probe velocities such that they all have the same radial velocity
vector2 *velocitySpectrum(vector2 *probePositions, Planet *sun, int numProbes);
// Returns an array structure of circular perturbations in xz plane read from file. Faster to read from array for use over multiple probes. Struct: data[NUMPOINTS][x, z]
vector2 *returnPerturbationsCircle(int NUMPERTURBED);
// Initializes all James Webb probes, their perturbed probes, and their h2 probes.
void initializeProbes(Probe *probes, vector2 *probePositions, vector2 *probeVels, vector2 *perturbedOffsets, int numProbes, int NUMPERTURBED, double RADIUS);

// RK4 for sun earth system with multiple James Webb Telescopes with perturbed probes, and h2 probes
void multiProbeRK4(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, double dt);
// Calculates and stores gravitational accelerations on earth, multiple probes, perturbations, and h2
void multiProbeCalcAccels(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels);
// Stores velocities in kVels
void multiProbeExtractVels(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthVel, vector2 *probeVels, vector2 **perturbedVels, vector2 **h2Vels);
// Moves Objects forward using kAccel and velocity
void multiProbeRK4Step(Planet *earth, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels, double dt);

// Lyapunov functions

// Handles lyapunov calculations at every timestep to calculate h1 and h2
void multiProbeLyapunovStuff(Probe *probes, double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double RADIUS, double initialArea);

// Extracts the lyapunov results from a multi probe simulation
void multiProbeExtractLyapunovs(double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double **lyapunovs, double t);





#endif /* HEADERS_MULTIPROBECHAOS_H_ */
