/*
 * circleGenerator.h
 *
 *  Created on: Feb 14, 2023
 *      Author: norms
 */

#ifndef HEADERS_CIRCLEGENERATOR_H_
#define HEADERS_CIRCLEGENERATOR_H_

#include <string>
#include "vector2.h"


void toFile(std::string fileName, double **data, int dataLen);
void writeCircleToFile(int numPoints, double radius);
void uniformCircle(int numPoints, double **arr, double r);


#endif /* HEADERS_CIRCLEGENERATOR_H_ */
