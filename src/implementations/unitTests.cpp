/*
 * unitTests.cpp
 *
 *  Created on: Feb 19, 2023
 *      Author: norms
 */

#include <iostream>
#include <iomanip>
#include <math.h>

#include "../headers/unitTests.h"


// Multi Probe Unit Test Functions

// Prints out the perturbed offsets from the unit circle generation
void perturbedOffsetTest(vector2 *perturbedOffsets, int NUMPERTURBED)
{
	std::cout << "Printing perturbed offsets" << std::endl;
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		std::cout << "Perturbed probe: " << i << std::endl;
		std::cout << "x: " << std::setprecision(10) << perturbedOffsets[i].x << " y: " << std::setprecision(10) << perturbedOffsets[i].y << std::endl << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Prints out the James Webb positions that populate the grid with specified width and height
void probePositionsTest(vector2 *probePositions, int numProbes, int width, int height)
{
	std::cout << "Printing probe positions" << std::endl;
	for (int i = 0; i < width * height; ++i)
	{
		std::cout << "Probe: " << i << std::endl;
		std::cout << "x: " << std::setprecision(10) << probePositions[i].x << " y: " << std::setprecision(10) << probePositions[i].y << std::endl << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Returns a vector2 array of probe velocities such that they all have the same radial velocity
vector2 *velocitySpectrumTest(vector2 *probePositions, Planet *sun, int numProbes)
{
	// Allocate space for probe velocities
	vector2 *probeVels = (vector2 *)malloc(numProbes * sizeof(vector2));

	// Loop through all probes
	std::cout << "Printing velocities for each probe" << std::endl;

	for (int i = 0; i < numProbes; ++i)
	{
		std::cout << "Probe " << i << std::endl;

		// Vector from the probe to the sun
		vector2 r = {sun->pos.x - probePositions[i].x, sun->pos.y - probePositions[i].y};

		// Unit vector pointing from probe to sun
		vector2 rUnit = unitVector(r);

		// Magnitude
		double dist = mag(r);

		std::cout << "Vector to sun is: x: " << std::setprecision(10) << r.x << " y: " << std::setprecision(10) << r.y << std::endl;
		std::cout << "Unit vector is: x: " << std::setprecision(10) << rUnit.x << " y: " << std::setprecision(10) << rUnit.y << std::endl;
		std::cout << "Distance to sun: " << std::setprecision(10) << dist << std::endl << std::endl;

		std::cout << "Calculating the equivalent radial velocity to the earth" << std::endl;

		// Radial speed of the earth: radialSpeed = v / earthDist
		double earthSpeed = 29787.7;
		double earthToSunDist = 149597870700;
		double radialSpeed = earthSpeed / earthToSunDist;

		// Find speed of probe: v = dist * radialSpeed
		double probeSpeed = radialSpeed * dist;

		// Perpendicular vector to unit vector in the counter clockwise direction
		vector2 rUnitOrbital = {-1 * rUnit.y, rUnit.x};

		std::cout << "Radial speed: " << std::setprecision(10) << radialSpeed << std::endl;
		std::cout << "Equivalent probe speed: " << std::setprecision(10) << probeSpeed << std::endl;
		std::cout << "Unit vector in orbital dir: x: " << std::setprecision(10) << rUnitOrbital.x << " y: " << std::setprecision(10) << rUnitOrbital.y << std::endl;

		// Calculate and store velocity vector
		probeVels[i] = {probeSpeed * rUnitOrbital.x, probeSpeed * rUnitOrbital.y};

		std::cout << "Calculated probe velocity: x: " << std::setprecision(10) << probeVels[i].x << " y: " << std::setprecision(10) << probeVels[i].y << std::endl << std::endl;
	}

	std::cout << std::endl << std::endl << "--------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
	// Return the probe velocity array
	return probeVels;
}

// Initializes all James Webb probes, their perturbed probes, and their h2 probes.
void initializeProbesTest(Probe *probes, vector2 *probePositions, vector2 *probeVels, vector2 *perturbedOffsets, int numProbes, int NUMPERTURBED, double RADIUS)
{
	std::cout << "Printing information of initialized probes, perturbed probes, and h2 probes" << std::endl << std::endl;
	for (int i = 0; i < numProbes; ++i)
	{
		// James Webb:
		// Position
		probes[i].pos = {probePositions[i].x, probePositions[i].y};

		// Velocity
		probes[i].vel = {probeVels[i].x, probeVels[i].y};

		std::cout << "James Webb Probe: " << i << std::endl;
		std::cout << "Pos: x: " << std::setprecision(10) << probes[i].pos.x << " y: " << std::setprecision(10) << probes[i].pos.y << std::endl;
		std::cout << "Vel: x: " << std::setprecision(10) << probes[i].vel.x << " y: " << std::setprecision(10) << probes[i].vel.y << std::endl << std::endl;
		std::cout << "---------------------------------" << std::endl << std::endl;
		// Perturbed Probes:
		// Allocate space for perturbed probes
		probes[i].perturbed = (Probe *)malloc(NUMPERTURBED * sizeof(Probe));

		std::cout << "Perturbed probes for probe " << i << std::endl;
		// Copy information over
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			std::cout << "Perturbed probe: " << j << std::endl;

			// Position: Add the perturbed offset to the James Webb position
			probes[i].perturbed[j].pos = {probes[i].pos.x + perturbedOffsets[j].x, probes[i].pos.y + perturbedOffsets[j].y};

			// Velocity:
			probes[i].perturbed[j].vel = {probes[i].vel.x, probes[i].vel.y};

			std::cout << "Pos: x: " << std::setprecision(10) << probes[i].perturbed[j].pos.x << " y: " << std::setprecision(10) << probes[i].perturbed[j].pos.y << std::endl;
			std::cout << "Vel: x: " << std::setprecision(10) << probes[i].perturbed[j].vel.x << " y: " << std::setprecision(10) << probes[i].perturbed[j].vel.y << std::endl << std::endl;
		}
		std::cout << "---------------------------------" << std::endl << std::endl;

		// h2 Orthonormal Probes:

		// Allocate space
		probes[i].h2 = (Probe *)malloc(2 * sizeof(Probe));

		std::cout << "h2 orthonormal probes" << std::endl;

		// Copy information over
		for (int j = 0; j < 2; ++j)
		{
			copyProbe(&probes[i].h2[j], probes[i]);
		}

		// Add orthonormal offsets
		probes[i].h2[0].pos.x += RADIUS;
		probes[i].h2[1].pos.y += RADIUS;

		std::cout << "Probe 0" << std::endl;
		std::cout << "Pos: x: " << std::setprecision(10) << probes[i].h2[0].pos.x << " y: " << std::setprecision(10) << probes[i].h2[0].pos.y << std::endl;
		std::cout << "Vel: x: " << std::setprecision(10) << probes[i].h2[0].vel.x << " y: " << std::setprecision(10) << probes[i].h2[0].vel.y << std::endl << std::endl;

		std::cout << "Probe 1" << std::endl;
		std::cout << "Pos: x: " << std::setprecision(10) << probes[i].h2[1].pos.x << " y: " << std::setprecision(10) << probes[i].h2[1].pos.y << std::endl;
		std::cout << "Vel: x: " << std::setprecision(10) << probes[i].h2[1].vel.x << " y: " << std::setprecision(10) << probes[i].h2[1].vel.y << std::endl << std::endl;
	}
	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Moves Objects forward using kAccel and velocity
void multiProbeRK4StepTest(Planet *earth, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels, double dt)
{
	std::cout << "Unit test for rk4 step function" << std::endl << std::endl;
	// Push velocities forward using k Accel, and then push positions forward using current velocity

	// Earth:
	// Velocities
	std::cout << "dt is: " << dt << std::endl;
	std::cout << "Earth" << std::endl;
	std::cout << "Old velocities: vx: " << std::setprecision(15) << earth->vel.x << " vy: " << std::setprecision(15) << earth->vel.y << std::endl;

	earth->vel.x += earthAccel->x * dt;
	earth->vel.y += earthAccel->y * dt;

	std::cout << "New velocities: vx: " << std::setprecision(15) << earth->vel.x << " vy: " << std::setprecision(15) << earth->vel.y << std::endl;
	std::cout << "Old positions: x: " << std::setprecision(15) << earth->pos.x << " y: " << std::setprecision(15) << earth->pos.y << std::endl;


	// Positions
	earth->pos.x += earth->vel.x * dt;
	earth->pos.y += earth->vel.y * dt;

	std::cout << "New positions: x: " << std::setprecision(15) << earth->pos.x << " y: " << std::setprecision(15) << earth->pos.y << std::endl << std::endl;

	// Probes:
	std::cout << "James Webb Probes" << std::endl << std::endl;
	for (int i = 0; i < numProbes; ++i)
	{
		// Velocities

		std::cout << "Probe: " << i << std::endl;
		std::cout << "Old velocities: vx: " << std::setprecision(15) << probes[i].vel.x << " vy: " << std::setprecision(15) << probes[i].vel.y << std::endl;

		probes[i].vel.x += probeAccels[i].x * dt;
		probes[i].vel.y += probeAccels[i].y * dt;

		std::cout << "New velocities: vx: " << std::setprecision(15) << probes[i].vel.x << " vy: " << std::setprecision(15) << probes[i].vel.y << std::endl;
		std::cout << "Old positions: x: " << std::setprecision(15) << probes[i].pos.x << " y: " << std::setprecision(15) << probes[i].pos.y << std::endl;

		// Positions
		probes[i].pos.x += probes[i].vel.x * dt;
		probes[i].pos.y += probes[i].vel.y * dt;

		std::cout << "New positions: x: " << std::setprecision(15) << probes[i].pos.x << " y: " << std::setprecision(15) << probes[i].pos.y << std::endl << std::endl;

		// Perturbed Probes:
		std::cout << "Perturbed Probes" << std::endl << std::endl;
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Velocities

			std::cout << "Perturbed Probe: " << j << std::endl;
			std::cout << "Old velocities: vx: " << std::setprecision(15) << probes[i].perturbed[j].vel.x << " vy: " << std::setprecision(15) << probes[i].perturbed[j].vel.y << std::endl;

			probes[i].perturbed[j].vel.x += perturbedAccels[i][j].x * dt;
			probes[i].perturbed[j].vel.y += perturbedAccels[i][j].y * dt;

			std::cout << "New velocities: vx: " << std::setprecision(15) << probes[i].perturbed[j].vel.x << " vy: " << std::setprecision(15) << probes[i].perturbed[j].vel.y << std::endl;
			std::cout << "Old positions: x: " << std::setprecision(15) << probes[i].perturbed[j].pos.x << " y: " << std::setprecision(15) << probes[i].perturbed[j].pos.y << std::endl;

			// Positions
			probes[i].perturbed[j].pos.x += probes[i].perturbed[j].vel.x * dt;
			probes[i].perturbed[j].pos.y += probes[i].perturbed[j].vel.y * dt;

			std::cout << "New positions: x: " << std::setprecision(15) << probes[i].perturbed[j].pos.x << " y: " << std::setprecision(15) << probes[i].perturbed[j].pos.y << std::endl << std::endl;
		}

		// h2 Probes:
		std::cout << "h2 Probes:" << std::endl << std::endl;
		for (int j = 0; j < 2; ++j)
		{
			// Velocities

			std::cout << "h2 Probe: " << j << std::endl;
			std::cout << "Old velocities: vx: " << std::setprecision(15) << probes[i].h2[j].vel.x << " vy: " << std::setprecision(15) << probes[i].h2[j].vel.y << std::endl;

			probes[i].h2[j].vel.x += h2Accels[i][j].x * dt;
			probes[i].h2[j].vel.y += h2Accels[i][j].y * dt;

			std::cout << "New velocities: vx: " << std::setprecision(15) << probes[i].h2[j].vel.x << " vy: " << std::setprecision(15) << probes[i].h2[j].vel.y << std::endl;
			std::cout << "Old positions: x: " << std::setprecision(15) << probes[i].h2[j].pos.x << " y: " << std::setprecision(15) << probes[i].h2[j].pos.y << std::endl;

			// Positions
			probes[i].h2[j].pos.x += probes[i].h2[j].vel.x * dt;
			probes[i].h2[j].pos.y += probes[i].h2[j].vel.y * dt;

			std::cout << "New positions: x: " << std::setprecision(15) << probes[i].h2[j].pos.x << " y: " << std::setprecision(15) << probes[i].h2[j].pos.y << std::endl << std::endl;
		}
	}
	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Stores velocities in kVels
void multiProbeExtractVelsTest(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthVel, vector2 *probeVels, vector2 **perturbedVels, vector2 **h2Vels)
{
	std::cout << "Unit test for velocity extraction" << std::endl << std::endl;
	// Earth
	earthVel->x = earth->vel.x;
	earthVel->y = earth->vel.y;

	std::cout << "Printing earth velocities" << std::endl;
	std::cout << "Original: vx: " << std::setprecision(10) << earth->vel.x << " vy: " << std::setprecision(10) << earth->vel.y << std::endl;
	std::cout << "Stored: vx: " << std::setprecision(10) << earthVel->x << " vy: " << std::setprecision(10) << earthVel->y << std::endl << std::endl;

	// Probes:
	std::cout << "James Webb Probes" << std::endl << std::endl;
	for (int i = 0; i < numProbes; ++i)
	{
		probeVels[i].x = probes[i].vel.x;
		probeVels[i].y = probes[i].vel.y;

		std::cout << "Printing Probe: " << i << " Velocities" << std::endl;
		std::cout << "Original: vx: " << std::setprecision(10) << probes[i].vel.x << " vy: " << std::setprecision(10) << probes[i].vel.y << std::endl;
		std::cout << "Stored: vx: " << std::setprecision(10) << probeVels[i].x << " vy: " << probeVels[i].y << std::endl << std::endl;

		// Perturbed Probes:
		std::cout << "Perturbed probes" << std::endl << std::endl;
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			perturbedVels[i][j].x = probes[i].perturbed[j].vel.x;
			perturbedVels[i][j].y = probes[i].perturbed[j].vel.y;

			std::cout << "Printing Perturbed Probe: " << j << " Velocities" << std::endl;
			std::cout << "Original: vx: " << std::setprecision(10) << probes[i].perturbed[j].vel.x << " vy: " << std::setprecision(10) << probes[i].perturbed[j].vel.y << std::endl;
			std::cout << "Stored: vx: " << std::setprecision(10) << perturbedVels[i][j].x << " vy: " << perturbedVels[i][j].y << std::endl << std::endl;
		}

		// h2 Probes:
		std::cout << "h2 probes" << std::endl << std::endl;
		for (int j = 0; j < 2; ++j)
		{
			h2Vels[i][j].x = probes[i].h2[j].vel.x;
			h2Vels[i][j].y = probes[i].h2[j].vel.y;

			std::cout << "Printing h2 Probe: " << j << " Velocities" << std::endl;
			std::cout << "Original: vx: " << std::setprecision(10) << probes[i].h2[j].vel.x << " vy: " << std::setprecision(10) << probes[i].h2[j].vel.y << std::endl;
			std::cout << "Stored: vx: " << std::setprecision(10) << h2Vels[i][j].x << " vy: " << h2Vels[i][j].y << std::endl << std::endl;
		}
	}
	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Calculates and stores gravitational accelerations on earth, multiple probes, perturbations, and h2
void multiProbeCalcAccelsTest(Planet *earth, Planet *sun, Probe *probes, int numProbes, int NUMPERTURBED, vector2 *earthAccel, vector2 *probeAccels, vector2 **perturbedAccels, vector2 **h2Accels)
{
	// Gravitational Constant
	double G = 6.6743e-11;

	// Sun on earth:
	vector2 r = {earth->pos.x - sun->pos.x, earth->pos.y - sun->pos.y}; 	// Vector pointing from sun to earth
	double magr = mag(r);													// Magnitude of that vector
	vector2 rUnit = unitVector(r);											// Direction of the vector
	double g = (-1 * G * sun->m) / (magr * magr);							// Newtons law of gravity to get force magnitude
	// Multiply by unit vector for direction
	earthAccel->x = g * rUnit.x;
	earthAccel->y = g * rUnit.y;
	std::cout << "Calc Accel Unit Test!" << std::endl << std::endl;

	std::cout << "Sun on earth" << std::endl;
	std::cout << "r: " << r.x << " " << r.y << std::endl;
	std::cout << "Magnitude of r: " << magr << std::endl;
	std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
	std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
	std::cout << "Acceleration: x: " << earthAccel->x << " y: " << earthAccel->y << std::endl << std::endl;

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

		std::cout << "Sun on Probe: " << i << std::endl;
		std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
		std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
		std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
		std::cout << "Acceleration: x: " << std::setprecision(10) << probeAccels[i].x << " y: " << std::setprecision(10) << probeAccels[i].y << std::endl << std::endl;

		// Earth on James Webb Probe
		r = {probes[i].pos.x - earth->pos.x, probes[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		// Add to accelerations from the Sun
		probeAccels[i].x += g * rUnit.x;
		probeAccels[i].y += g * rUnit.y;

		std::cout << "Earth on Probe: " << i << std::endl;
		std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
		std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
		std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
		std::cout << "Acceleration from earth: x: " << std::setprecision(10) << g * rUnit.x << " y: " << std::setprecision(10) << g * rUnit.y << std::endl;
		std::cout << "Net Acceleration: x: " << std::setprecision(10) << probeAccels[i].x << " y: " << std::setprecision(10) << probeAccels[i].y << std::endl << std::endl;

		std::cout << "Perturbed Probes:" << std::endl << std::endl;

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

			std::cout << "Sun on perturbed probe: " << j << std::endl;
			std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
			std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
			std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
			std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
			std::cout << "Acceleration: x: " << std::setprecision(10) << perturbedAccels[i][j].x << " y: " << std::setprecision(10) << perturbedAccels[i][j].y << std::endl << std::endl;

			// Earth on perturbation
			r = {probes[i].perturbed[j].pos.x - earth->pos.x, probes[i].perturbed[j].pos.y - earth->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * earth->m) / (magr * magr);
			perturbedAccels[i][j].x += g * rUnit.x;
			perturbedAccels[i][j].y += g * rUnit.y;

			std::cout << "Earth on perturbed probe: " << j << std::endl;
			std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
			std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
			std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
			std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
			std::cout << "Acceleration from earth: x: " << std::setprecision(10) << g * rUnit.x << " y: " << std::setprecision(10) << g * rUnit.y << std::endl;
			std::cout << "Net Acceleration: x: " << std::setprecision(10) << perturbedAccels[i][j].x << " y: " << std::setprecision(10) << perturbedAccels[i][j].y << std::endl << std::endl;
		}

		std::cout << "h2 Probes: " << std::endl << std::endl;
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

			std::cout << "Sun on h2 probe: " << j << std::endl;
			std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
			std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
			std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
			std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
			std::cout << "Acceleration: x: " << std::setprecision(10) << h2Accels[i][j].x << " y: " << std::setprecision(10) << h2Accels[i][j].y << std::endl << std::endl;

			// Earth on h2 probe
			r = {probes[i].h2[j].pos.x - earth->pos.x, probes[i].h2[j].pos.y - earth->pos.y};
			magr = mag(r);
			rUnit = unitVector(r);
			g = (-1 * G * earth->m) / (magr * magr);
			h2Accels[i][j].x += g * rUnit.x;
			h2Accels[i][j].y += g * rUnit.y;

			std::cout << "Earth on h2 probe: " << j << std::endl;
			std::cout << "r: " << std::setprecision(10) << r.x << " " << std::setprecision(10) << r.y << std::endl;
			std::cout << "Magnitude of r: " << std::setprecision(10) << magr << std::endl;
			std::cout << "Unit vector for r: " << std::setprecision(10) << rUnit.x << " " << std::setprecision(10) << rUnit.y << std::endl;
			std::cout << "Magnitude of gravitational acceleration: " << std::setprecision(10) << g << std::endl;
			std::cout << "Acceleration from earth: x: " << std::setprecision(10) << g * rUnit.x << " y: " << std::setprecision(10) << g * rUnit.y << std::endl;
			std::cout << "Net Acceleration: x: " << std::setprecision(10) << h2Accels[i][j].x << " y: " << std::setprecision(10) << h2Accels[i][j].y << std::endl << std::endl;
		}
	}
	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

// Extracts the lyapunov results from a multi probe simulation
void multiProbeExtractLyapunovsTest(double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double **lyapunovs, double t)
{
	std::cout << "Unit test for lyapunov extraction function for multi probe simulation" << std::endl << std::endl;
	// Loop through probes
	for (int i = 0; i < numProbes; ++i)
	{
		std::cout << "Probe: " << i << std::endl << std::endl;
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

		std::cout << "Maximum LfSum is: " << max << std::endl;

		// The maximum lyapunov exponent is the maximum divided by the total elapsed time
		lyapunovs[i][0] = max / t;

		std::cout << "h1 lyapunov exponent is: " << lyapunovs[i][0] << std::endl;

		// Calculate and store the second lyapunov exponent
		lyapunovs[i][1] = (areaSums[i] / t) - lyapunovs[i][0];

		std::cout << "Area sum is: " << areaSums[i] << std::endl;

		std::cout << "h2 lyapunov exponent is: " << lyapunovs[i][1] << std::endl;
	}
}

// Handles lyapunov calculations at every timestep to calculate h1 and h2
void multiProbeLyapunovStuffTest(Probe *probes, double **LfSums, double *areaSums, int numProbes, int NUMPERTURBED, double RADIUS, double initialArea)
{
	// h1: Max lyapunov
	std::cout << "Unit test for timestep calculations for lyapunov exponents" << std::endl << std::endl;
	std::cout << "For h1" << std::endl << std::endl;
	// Loop through all probes and perturbed Probes
	for (int i = 0; i < numProbes; ++i)
	{
		std::cout << std::endl << std::endl << "Probe: " << i << std::endl << std::endl;

		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// Calculate Lf of each perturbation, and renormalize the position vector
			// Vector pointing from probe to perturbed probe
			vector2 v = {probes[i].perturbed[j].pos.x - probes[i].pos.x, probes[i].perturbed[j].pos.y - probes[i].pos.y};

			// Find the distance of the vector and scale it relative to the inital distance
			double dist = mag(v) / RADIUS;

			// Take the natural log
			double logDist = log(dist);

			std::cout << std::endl << "Perturbed probe: " << j << std::endl;
			std::cout << "Distance from james webb is: " << std::setprecision(15) << dist * RADIUS << std::endl;
			std::cout << "Scaled distance: " << std::setprecision(15) << dist << std::endl;


			// Add to the running sum for this perturbed probe

			std::cout << "Lfsum before was: " << std::setprecision(15) << LfSums[i][j] << std::endl;

			LfSums[i][j] += logDist;

			std::cout << "Lfsum after is: " << std::setprecision(15) << LfSums[i][j] << std::endl;

			// Renormalize the perturbed probe preserving the direction vector
			vector2 uv = unitVector(v);

			std::cout << "Vector from james webb to probe" << std::endl;
			std::cout << "x: " << std::setprecision(15) << v.x << " y: " << std::setprecision(15) << v.y << std::endl;


			// Scale unit vector by the initial distance
			uv.x *= RADIUS;
			uv.y *= RADIUS;

			std::cout << "Scaled unit vector from james webb to probe" << std::endl;
			std::cout << "x: " << std::setprecision(15) << uv.x << " y: " << std::setprecision(15) << uv.y << std::endl;

			// Renormalize the perturbed probe

			std::cout << "Perturbed original pos: x: " << std::setprecision(15) << probes[i].perturbed[j].pos.x << " y: " << std::setprecision(15) << probes[i].perturbed[j].pos.y << std::endl;
			std::cout << "Probe's position is: x: " << std::setprecision(15) << probes[i].pos.x << " y: " << std::setprecision(15) << probes[i].pos.y << std::endl;
			probes[i].perturbed[j].pos = {probes[i].pos.x + uv.x, probes[i].pos.y + uv.y};

			std::cout << "Perturbed new pos: x: " << std::setprecision(15) << probes[i].perturbed[j].pos.x << " y: " << std::setprecision(15) << probes[i].perturbed[j].pos.y << std::endl;
		}
		std::cout << std::endl;

		// h2 lyapunov calculation:

		std::cout << "For h2" << std::endl << std::endl;

		// Two vectors pointing from probe to each h2 probe
		vector2 v1 = {probes[i].h2[0].pos.x - probes[i].pos.x, probes[i].h2[0].pos.y - probes[i].pos.y};
		vector2 v2 = {probes[i].h2[1].pos.x - probes[i].pos.x, probes[i].h2[1].pos.y - probes[i].pos.y};

		std::cout << "Vector from probe to first h2 probe" << std::endl;
		std::cout << "x: " << v1.x << std::setprecision(15) << " y: " << std::setprecision(15) << v1.y << std::endl;
		std::cout << "Vector from probe to second h2 probe" << std::endl;
		std::cout << "x: " << std::setprecision(15) << v2.x << " y: " << std::setprecision(15) << v2.y << std::endl << std::endl;

		// Calculate area of parallelogram formed by two vectors: A = ||a|| ||b|| sin(ab)
		double area = mag(v1) * mag(v2) * sin(angle(v1, v2));

		// Scale the area relative to the initial area
		double scaledArea = area / initialArea;

		std::cout << "Area of the parallelogram is: " << std::setprecision(15) << area << std::endl;
		std::cout << "Divided by initial area: " << std::setprecision(15) << scaledArea << std::endl;

		// Take natural log
		double logArea = log(scaledArea);

		std::cout << "Log area of the parallelogram is: " << std::setprecision(15) << logArea << std::endl;

		// Add to the areaSum

		std::cout << "Area sum before was: " << std::setprecision(15) << areaSums[i] << std::endl;

		areaSums[i] += logArea;

		std::cout << "Area sum after is: " << std::setprecision(15) << areaSums[i] << std::endl << std::endl;

		// Renormalize the vectors using Gram Schmidt Algorithm to preserve direction:
		// First orthonormal basis vector

		std::cout << "Gram schmidt orthogonalization algorithm" << std::endl << std::endl;

		std::cout << "Original vectors" << std::endl;
		std::cout << "V1: x: " << std::setprecision(15) << v1.x << " y: " << std::setprecision(15) << v1.y << std::endl;
		std::cout << "V2: x: " << std::setprecision(15) << v2.x << " y: " << std::setprecision(15) << v2.y << std::endl << std::endl;

		vector2 e1 = unitVector(v1);

		// Dot product between v2 and e1
		double v2Dote1 = dotProd(v2, e1);

		// Calculate the orthogonal basis vector
		vector2 y2 = {v2.x - (v2Dote1) * e1.x, v2.y - (v2Dote1) * e1.y};

		// Take unit vector to get second orthonormal basis vector
		vector2 e2 = unitVector(y2);


		std::cout << "Calculated orthonormal vectors" << std::endl;
		std::cout << "e1: " << std::setprecision(15) << e1.x << " " << e1.y << std::endl;
		std::cout << "e2: " << std::setprecision(15) << e2.x << " " << e2.y << std::endl << std::endl;

		// Scale the basis vectors by the initial distance
		e1.x *= RADIUS;
		e1.y *= RADIUS;
		e2.x *= RADIUS;
		e2.y *= RADIUS;

		std::cout << "Scaled orthonormal vectors" << std::endl;
		std::cout << "e1: " << std::setprecision(15) << e1.x << " " << e1.y << std::endl;
		std::cout << "e2: " << std::setprecision(15) << e2.x << " " << e2.y << std::endl << std::endl;

		std::cout << "Assigning the new vectors" << std::endl;
		std::cout << "Old positions" << std::endl;
		std::cout << "h2 Probe 1: x: " << std::setprecision(15) << probes[i].h2[0].pos.x << " y: " << std::setprecision(15) << probes[i].h2[0].pos.y << std::endl;
		std::cout << "h2 Probe 2: x: " << std::setprecision(15) << probes[i].h2[1].pos.x << " y: " << std::setprecision(15) << probes[i].h2[1].pos.y << std::endl;

		// Renormalize the h2 probes
		probes[i].h2[0].pos = {probes[i].pos.x + e1.x, probes[i].pos.y + e1.y};
		probes[i].h2[1].pos = {probes[i].pos.x + e2.x, probes[i].pos.y + e2.y};

		std::cout << "New positions" << std::endl;
		std::cout << "h2 Probe 1: x: " << std::setprecision(15) << probes[i].h2[0].pos.x << " y: " << std::setprecision(15) << probes[i].h2[0].pos.y << std::endl;
		std::cout << "h2 Probe 2: x: " << std::setprecision(15) << probes[i].h2[1].pos.x << " y: " << std::setprecision(15) << probes[i].h2[1].pos.y << std::endl;
	}
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Single Probe Unit Test Functions
// Lyapunov unit tests

// Gram schmidt to re-orthonormalize unit vectors
void gramSchmidtRenormalizationh2Test(Probe *jw, int NUMPERTURBED, vector2 *v1, vector2 *v2, double RADIUS)
{
	// Calculate orthonormal vectors e1 and e2
	vector2 e1 = unitVector(*v1);											// First orthogonal normalized basis vector

	std::cout << "Gram Schmidt Renormalization" << std::endl;
	std::cout << "Original vectors" << std::endl;
	std::cout << "V1: x: " << std::setprecision(15) << v1->x << " y: " << std::setprecision(15) << v1->y << std::endl;
	std::cout << "V2: x: " << std::setprecision(15) << v2->x << " y: " << std::setprecision(15) << v2->y << std::endl << std::endl;

	double v2Dote1 = dotProd(*v2, e1);										// Dot prodoct of second vector with orthobasis vector
	vector2 y2 = {v2->x - (v2Dote1) * e1.x, v2->y - (v2Dote1) * e1.y};		// Get the second orthogonal basis vector
	vector2 e2 = unitVector(y2);											// Normalize to get second orthonormal vector

	// Scale orthormal vectors by the radius
	e1.x *= RADIUS;
	e1.y *= RADIUS;
	e2.x *= RADIUS;
	e2.y *= RADIUS;


	std::cout << "New scaled orthonormal vectors" << std::endl;
	std::cout << "e1: " << std::setprecision(15) << e1.x << " " << e1.y << std::endl;
	std::cout << "e2: " << std::setprecision(15) << e2.x << " " << e2.y << std::endl << std::endl;

	// Assign new vectors
	std::cout << "Assigning the new vectors" << std::endl;
	std::cout << "Old positions" << std::endl;
	std::cout << "h2 Probe 1: x: " << std::setprecision(15) << jw->h2[0].pos.x << " y: " << std::setprecision(15) << jw->h2[0].pos.y << std::endl;
	std::cout << "h2 Probe 2: x: " << std::setprecision(15) << jw->h2[1].pos.x << " y: " << std::setprecision(15) << jw->h2[1].pos.y << std::endl;

	jw->h2[0].pos.x = jw->pos.x + e1.x; jw->h2[0].pos.y = jw->pos.y + e1.y;
	jw->h2[1].pos.x = jw->pos.x + e2.x; jw->h2[1].pos.y = jw->pos.y + e2.y;

	std::cout << "New positions" << std::endl;
	std::cout << "h2 Probe 1: x: " << std::setprecision(15) << jw->h2[0].pos.x << " y: " << std::setprecision(15) << jw->h2[0].pos.y << std::endl;
	std::cout << "h2 Probe 2: x: " << std::setprecision(15) << jw->h2[1].pos.x << " y: " << std::setprecision(15) << jw->h2[1].pos.y << std::endl << std::endl;
}

// Extracts lyapunov exponents for single probe with perturbations
void extractLyapunovsTest(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t)
{
	std::cout << "Unit test for lyapunov extraction function for single probe simulation" << std::endl << std::endl;
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

	std::cout << "Maximum LfSum is: " << max << std::endl;

	// The maximum lyapunov exponent is the maximum divided by the total elapsed time
	lyapunovs[0] = max / t;

	std::cout << "h1 lyapunov exponent is: " << lyapunovs[0] << std::endl;

	// Calculate and store the second lyapunov exponent
	lyapunovs[1] = (areaSum / t) - lyapunovs[0];

	std::cout << "Area sum is: " << areaSum << std::endl;
	std::cout << "h2 lyapunov exponent is: " << lyapunovs[1] << std::endl;
}

void lyapunovChaosStuffTest(Probe *jw, int NUMPERTURBED, double *LfSums, double initialDist, double *areaSum, double initialArea, double RADIUS)
{
	// h1: Max lyapunov
	std::cout << "For h1" << std::endl;
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Calculate Lf of each perturbation, and renormalize the position vector
		vector2 v = {jw->perturbed[i].pos.x - jw->pos.x, jw->perturbed[i].pos.y - jw->pos.y};

		// Take the current distance
		double dist = mag(v);

		// Divide current distance by original distance
		double scaledDist = dist / RADIUS;

		// Take the log
		double logDist = log(scaledDist);

		std::cout << "Perturbed probe: " << i << std::endl;
		std::cout << "Distance from james webb is: " << std::setprecision(15) << dist << " Log is: " << std::setprecision(15) << logDist << std::endl;
		std::cout << "Divided by initial distance: " << std::setprecision(15) << scaledDist << std::endl;

		// Add to the running sum
		std::cout << "Running sum before was: " << std::setprecision(15) << LfSums[i] << std::endl;
		LfSums[i] += logDist;

		std::cout << "Running sum after is: " << std::setprecision(15) << LfSums[i] << std::endl;

		// Renormalize the perturbed probe preserving direction to original distance
		std::cout << "Vector from james webb to probe" << std::endl;
		std::cout << "x: " << std::setprecision(15) << v.x << " y: " << std::setprecision(15) << v.y << std::endl;

		vector2 uv = unitVector(v);
		uv.y *= RADIUS;
		uv.x *= RADIUS;


		std::cout << "Scaled unit vector from james webb to probe" << std::endl;
		std::cout << "x: " << std::setprecision(15) << uv.x << " y: " << std::setprecision(15) << uv.y << std::endl;

		std::cout << "Perturbed original pos: x: " << std::setprecision(15) << jw->perturbed[i].pos.x << " y: " << std::setprecision(15) << jw->perturbed[i].pos.y << std::endl;
		jw->perturbed[i].pos.x = jw->pos.x + uv.x; jw->perturbed[i].pos.y = jw->pos.y + uv.y;

		std::cout << "After normalization: x:" << std::setprecision(15) << jw->perturbed[i].pos.x << " y: " << std::setprecision(15) << jw->perturbed[i].pos.y << std::endl << std::endl;
	}
	std::cout << std::endl << std::endl;

	// h2 Lyapunov
	std::cout << "h2 Lyapunov" << std::endl;
	vector2 v1 = {jw->h2[0].pos.x - jw->pos.x, jw->h2[0].pos.y - jw->pos.y};
	vector2 v2 = {jw->h2[1].pos.x - jw->pos.x, jw->h2[1].pos.y - jw->pos.y};

	std::cout << "Printing vector for first h2 probe" << std::endl;
	std::cout << "x: " << v1.x << std::setprecision(15) << " y: " << std::setprecision(15) << v1.y << std::endl;
	std::cout << "Printing vector for second h2 probe" << std::endl;
	std::cout << "x: " << std::setprecision(15) << v2.x << " y: " << std::setprecision(15) << v2.y << std::endl;

	// Calculate area of parallelogram: A = ||a|| ||b|| sin(ab)
	double area = (mag(v1) * mag(v2) * sin(angle(v1, v2)));

	// Divide area by initial area
	double scaledArea = area / initialArea;

	std::cout << "Area of the parallelogram is: " << std::setprecision(15) << area << std::endl;
	std::cout << "Divided by initial area: " << std::setprecision(15) << scaledArea << std::endl;

	// Take the log
	double logArea = log(scaledArea);

	std::cout << "Log area of the parallelogram is: " << std::setprecision(15) << logArea << std::endl;
	// Add to the running sum
	std::cout << "Area sum before was: " << std::setprecision(15) << *areaSum << std::endl;
	(*areaSum) += logArea;
	std::cout << "Area sum after is: " << std::setprecision(15) << *areaSum << std::endl;
//	std::cout << std::setprecision(15) << area << std::endl;
//	std::cout << "magv1: " << mag(v1) << " magv2: " << mag(v2) << std::endl;

	// Reorthonormalize h2 vectors
	gramSchmidtRenormalizationh2Test(jw, NUMPERTURBED, &v1, &v2, RADIUS);

	std::cout << std::endl << std::endl << "------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}
// RK4 unit tests
void rk4StepTest(Planet *earth, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels, double dt)
{
	// Push velocities forward using k Accel, and then push positions forward using current velocity

	// Earth:
	// Velocities
	std::cout << "dt is: " << dt << std::endl;
	std::cout << "Earth" << std::endl;
	std::cout << "Old velocities: vx: " << std::setprecision(15) << earth->vel.x << " vy: " << std::setprecision(15) << earth->vel.y << std::endl;

	earth->vel.x += earthAccel->x * dt;
	earth->vel.y += earthAccel->y * dt;

	std::cout << "New velocities: vx: " << std::setprecision(15) << earth->vel.x << " vy: " << std::setprecision(15) << earth->vel.y << std::endl;
	std::cout << "Old positions: x: " << std::setprecision(15) << earth->pos.x << " y: " << std::setprecision(15) << earth->pos.y << std::endl;

	// Positions
	earth->pos.x += earth->vel.x * dt;
	earth->pos.y += earth->vel.y * dt;

	std::cout << "New positions: x: " << std::setprecision(15) << earth->pos.x << " y: " << std::setprecision(15) << earth->pos.y << std::endl;


	// James Webb:
	// Velocities

	std::cout << "James Webb" << std::endl;
	std::cout << "Old velocities: vx: " << std::setprecision(15) << jw->vel.x << " vy: " << std::setprecision(15) << jw->vel.y << std::endl;

	jw->vel.x += jwAccel->x * dt;
	jw->vel.y += jwAccel->y * dt;

	std::cout << "New velocities: vx: " << std::setprecision(15) << jw->vel.x << " vy: " << std::setprecision(15) << jw->vel.y << std::endl;
	std::cout << "Old positions: x: " << std::setprecision(15) << jw->pos.x << " y: " << std::setprecision(15) << jw->pos.y << std::endl;

	// Positions
	jw->pos.x += jw->vel.x * dt;
	jw->pos.y += jw->vel.y * dt;

	std::cout << "New positions: x: " << std::setprecision(15) << jw->pos.x << " y: " << std::setprecision(15) << jw->pos.y << std::endl;

	// Perturbations:
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Velocities
		std::cout << "Perturbed probe: " << i << std::endl;
		std::cout << "Old velocities: vx: " << std::setprecision(15) << jw->perturbed[i].vel.x << " vy: " << std::setprecision(15) << jw->perturbed[i].vel.y << std::endl;

		jw->perturbed[i].vel.x += perturbedAccels[i].x * dt;
		jw->perturbed[i].vel.y += perturbedAccels[i].y * dt;

		std::cout << "New velocities: vx: " << std::setprecision(15) << jw->perturbed[i].vel.x << " vy: " << std::setprecision(15) << jw->perturbed[i].vel.y << std::endl;
		std::cout << "Old positions: x: " << std::setprecision(15) << jw->perturbed[i].pos.x << " y: " << std::setprecision(15) << jw->perturbed[i].pos.y << std::endl;

		// Positions
		jw->perturbed[i].pos.x += jw->perturbed[i].vel.x * dt;
		jw->perturbed[i].pos.y += jw->perturbed[i].vel.y * dt;

		std::cout << "New positions: x: " << std::setprecision(15) << jw->perturbed[i].pos.x << " y: " << std::setprecision(15) << jw->perturbed[i].pos.y << std::endl;
	}

	// h2:
	for (int i = 0; i < 2; ++i)
	{
		// Velocities

		std::cout << "h2 probe: " << i << std::endl;
		std::cout << "Old velocities: vx: " << std::setprecision(15) << jw->h2[i].vel.x << " vy: " << std::setprecision(15) << jw->h2[i].vel.y << std::endl;

		jw->h2[i].vel.x += h2Accels[i].x * dt;
		jw->h2[i].vel.y += h2Accels[i].y * dt;

		std::cout << "New velocities: vx: " << std::setprecision(15) << jw->h2[i].vel.x << " vy: " << std::setprecision(15) << jw->h2[i].vel.y << std::endl;
		std::cout << "Old positions: x: " << std::setprecision(15) << jw->h2[i].pos.x << " y: " << std::setprecision(15) << jw->h2[i].pos.y << std::endl;

		// Positions
		jw->h2[i].pos.x += jw->h2[i].vel.x * dt;
		jw->h2[i].pos.y += jw->h2[i].vel.y * dt;

		std::cout << "New positions: x: " << std::setprecision(15) << jw->h2[i].pos.x << " y: " << std::setprecision(15) << jw->h2[i].pos.y << std::endl;
	}
	std::cout << std::endl << std::endl << std::endl << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}
void extractVelTest(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthVel, vector2 *jwVel, vector2 *perturbedVels, vector2 *h2Vels)
{
	// Earth
	earthVel->x = earth->vel.x;
	earthVel->y = earth->vel.y;

	// James Webb
	jwVel->x = jw->vel.x;
	jwVel->y = jw->vel.y;

	std::cout << "Printing earth velocities" << std::endl;
	std::cout << "Original: vx: " << earth->vel.x << " vy: " << earth->vel.y << std::endl;
	std::cout << "Stored: vx: " << earthVel->x << " vy: " << earthVel->y << std::endl << std::endl;

	std::cout << "Printing james webb velocities" << std::endl;
	std::cout << "Original: vx: " << jw->vel.x << " vy: " << jw->vel.y << std::endl;
	std::cout << "Stored: vx: " << jwVel->x << " vy: " << jwVel->y << std::endl << std::endl;

	// Perturbed
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		perturbedVels[i].x = jw->perturbed[i].vel.x;
		perturbedVels[i].y = jw->perturbed[i].vel.y;

		std::cout << "Printing perturbed " << i << " Velocities" << std::endl;
		std::cout << "Original: vx: " << jw->perturbed[i].vel.x << " vy: " << jw->perturbed[i].vel.y << std::endl;
		std::cout << "Stored: vx: " << perturbedVels[i].x << " vy: " << perturbedVels[i].y << std::endl << std::endl;
	}

	// h2
	for (int i = 0; i < 2; ++i)
	{
		h2Vels[i].x = jw->h2[i].vel.x;
		h2Vels[i].y = jw->h2[i].vel.y;

		std::cout << "Printing h2 probe " << i << " Velocities" << std::endl;
		std::cout << "Original: vx: " << jw->h2[i].vel.x << " vy: " << jw->h2[i].vel.y << std::endl;
		std::cout << "Stored: vx: " << h2Vels[i].x << " vy: " << h2Vels[i].y << std::endl << std::endl;
	}
	std::cout << std::endl << std::endl << std::endl << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}

void calcAccelTest(Planet *earth, Planet *sun, Probe *jw, int NUMPERTURBED, vector2 *earthAccel, vector2 *jwAccel, vector2 *perturbedAccels, vector2 *h2Accels)
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

	std::cout << "Sun on earth" << std::endl;
	std::cout << "r: " << r.x << " " << r.y << std::endl;
	std::cout << "Magnitude of r: " << magr << std::endl;
	std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
	std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
	std::cout << "Acceleration: x: " << earthAccel->x << " y: " << earthAccel->y << std::endl << std::endl;

	// James Webb Telescope: Add up accelerations from Sun and Earth: Same process

	// Sun on James Webb
	r = {jw->pos.x - sun->pos.x, jw->pos.y - sun->pos.y};
	magr = mag(r);
	rUnit = unitVector(r);
	g = (-1 * G * sun->m) / (magr * magr);
	jwAccel->x = g * rUnit.x;
	jwAccel->y = g * rUnit.y;

	std::cout << "Sun on James Webb" << std::endl;
	std::cout << "r: " << r.x << " " << r.y << std::endl;
	std::cout << "Magnitude of r: " << magr << std::endl;
	std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
	std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
	std::cout << "Acceleration: x: " << jwAccel->x << " y: " << jwAccel->y << std::endl << std::endl;

	// Earth on James Webb
	r = {jw->pos.x - earth->pos.x, jw->pos.y - earth->pos.y};
	magr = mag(r);
	rUnit = unitVector(r);
	g = (-1 * G * earth->m) / (magr * magr);
	// Add to accelerations from the Sun
	jwAccel->x += g * rUnit.x;
	jwAccel->y += g * rUnit.y;

	std::cout << "Earth on James Webb" << std::endl;
	std::cout << "r: " << r.x << " " << r.y << std::endl;
	std::cout << "Magnitude of r: " << magr << std::endl;
	std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
	std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
	std::cout << "Acceleration from earth: x: " << g * rUnit.x << " y: " << g * rUnit.y << std::endl;
	std::cout << "Net Acceleration: x: " << jwAccel->x << " y: " << jwAccel->y << std::endl << std::endl;

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

		std::cout << "Sun on Perturbation: " << i << std::endl;
		std::cout << "r: " << r.x << " " << r.y << std::endl;
		std::cout << "Magnitude of r: " << magr << std::endl;
		std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
		std::cout << "Acceleration: x: " << perturbedAccels[i].x << " y: " << perturbedAccels[i].y << std::endl << std::endl;

		// Earth on perturbation
		r = {jw->perturbed[i].pos.x - earth->pos.x, jw->perturbed[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		perturbedAccels[i].x += g * rUnit.x;
		perturbedAccels[i].y += g * rUnit.y;

		std::cout << "Earth on Perturbation: " << i << std::endl;
		std::cout << "r: " << r.x << " " << r.y << std::endl;
		std::cout << "Magnitude of r: " << magr << std::endl;
		std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
		std::cout << "Acceleration from earth: x: " << g * rUnit.x << " y: " << g * rUnit.y << std::endl;
		std::cout << "Net Acceleration: x: " << perturbedAccels[i].x << " y: " << perturbedAccels[i].y << std::endl << std::endl;
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


		std::cout << "Sun on h2 Probe: " << i << std::endl;
		std::cout << "r: " << r.x << " " << r.y << std::endl;
		std::cout << "Magnitude of r: " << magr << std::endl;
		std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
		std::cout << "Acceleration: x: " << h2Accels[i].x << " y: " << h2Accels[i].y << std::endl << std::endl;

		// Earth on h2 probe
		r = {jw->h2[i].pos.x - earth->pos.x, jw->h2[i].pos.y - earth->pos.y};
		magr = mag(r);
		rUnit = unitVector(r);
		g = (-1 * G * earth->m) / (magr * magr);
		h2Accels[i].x += g * rUnit.x;
		h2Accels[i].y += g * rUnit.y;

		std::cout << "Earth on h2 Probe: " << i << std::endl;
		std::cout << "r: " << r.x << " " << r.y << std::endl;
		std::cout << "Magnitude of r: " << magr << std::endl;
		std::cout << "Unit vector for r: " << rUnit.x << " " << rUnit.y << std::endl;
		std::cout << "Magnitude of gravitational acceleration: " << g << std::endl;
		std::cout << "Acceleration from earth: x: " << g * rUnit.x << " y: " << g * rUnit.y << std::endl;
		std::cout << "Net Acceleration: x: " << h2Accels[i].x << " y: " << h2Accels[i].y << std::endl << std::endl;
	}

	std::cout << std::endl << std::endl << std::endl << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}

