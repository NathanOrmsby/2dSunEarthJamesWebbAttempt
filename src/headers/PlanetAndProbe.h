/*
 * planetsAndProbe.h
 *
 *  Created on: Feb 11, 2023
 *      Author: norms
 */

#ifndef PLANETSANDPROBE_H_
#define PLANETSANDPROBE_H_

#include "vector2.h"

// Planets apply gravity
typedef struct
{
	vector2 pos; vector2 vel;
	double m;
} Planet;

// James Webb probe does not apply gravity. Contains an pointer for perturbation probes, and a pointer for probes used to calculate h2.
typedef struct Probe
{
	vector2 pos; vector2 vel;
	Probe *perturbed;
	Probe *h2;
} Probe;

// Utility functions for 2d Planets and Probes

// Planets
// Copies all data from i to o
void copyPlanet(Planet *o, Planet i);
// Probes
// Copies position and velocity data from one probe to another, doesn't worry about pointers
void copyProbe(Probe *o, Probe i);



#endif /* PLANETSANDPROBE_H_ */
