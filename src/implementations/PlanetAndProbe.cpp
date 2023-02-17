/*
 * PlanetAndProbe.cpp
 *
 *  Created on: Feb 14, 2023
 *      Author: norms
 */


#include "../headers/PlanetAndProbe.h"
#include "../headers/vector2.h"

// Utility functions for 2d Planets and Probes

// Planets
// Copies all data from i to o
void copyPlanet(Planet *o, Planet i)
{
	o->pos.x = i.pos.x; o->pos.y = i.pos.y; o->vel.x = i.vel.x; o->vel.y = i.vel.y; o->m = i.m;
}

// Probes
// Copies position and velocity data from one probe to another, doesn't worry about pointers
void copyProbe(Probe *o, Probe i)
{
	o->pos.x = i.pos.x; o->pos.y = i.pos.y; o->vel.x = i.vel.x; o->vel.y = i.vel.y;
}

