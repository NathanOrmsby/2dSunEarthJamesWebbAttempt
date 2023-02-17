/*
 * vector2.h
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_VECTOR2_H_
#define HEADERS_VECTOR2_H_


typedef struct
{
	double x;
	double y;
} vector2;

// Returns dot product of two vector2s
double dotProd(vector2 v1, vector2 v2);
// Returns magnitude of vector2
double mag(vector2 v);
// Returns unit vector of a vector2
vector2 unitVector(vector2 v);
// Returns angle between two vector2s
double angle(vector2 v1, vector2 v2);
// Returns angle between unit vectors of v1 and v2, good if worried about overflow
double unitAngle(vector2 v1, vector2 v2);



#endif /* HEADERS_VECTOR2_H_ */
