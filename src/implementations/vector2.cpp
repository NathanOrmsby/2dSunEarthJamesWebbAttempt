/*
 * vector2.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

// LINUX
#include "../headers/vector2.h"

// WINDOWS
//#include "..\headers\vector2.h"

#include "math.h"

double dotProd(vector2 v1, vector2 v2)
{
	return (v1.x * v2.x) + (v1.y * v2.y);
}

double mag(vector2 v)
{
	return sqrt((v.x * v.x) + (v.y * v.y));
}

// Returns unit vector of a vector2
vector2 unitVector(vector2 v)
{
	double vmag = mag(v);
	return (vector2){v.x / vmag, v.y / vmag};
}

// Returns angle between two vector2s
double angle(vector2 v1, vector2 v2)
{
	return acos(dotProd(v1, v2) / (mag(v1) * (mag(v2))));
}

// Returns angle between unit vectors of v1 and v2, good if worried about overflow
double unitAngle(vector2 v1, vector2 v2)
{
	vector2 uv1 = unitVector(v1); vector2 uv2 = unitVector(v2);
	return acos(dotProd(uv1, uv2));
}



