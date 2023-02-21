/*
 * toFile.h
 *
 *  Created on: Feb 19, 2023
 *      Author: norms
 */

#ifndef HEADERS_TOFILE_H_
#define HEADERS_TOFILE_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>

#include "vector2.h"

// Writes h2 probe evolution relative to james webb
void h2Evolution(vector2 *jwData, vector2 **h2Data, int dataLen);

// Writes evolution of circle of probes into ellipse
void perturbedEvolution(vector2 *jwData, vector2 **perturbedData, int dataLen, int NUMPERTURBED);
// Writes position data of perturbed probes relative to the james webb to csv. This function only works if there are four perturbed probes
void perturbedSeparation(vector2 *jwData, vector2 **perturbedData, int dataLen, int NUMPERTURBED);

// Writes position data to csv file for earth and james webb
void earthJWdataToFile(vector2 *earthData, vector2 *jwData, vector2 **perturbedData, vector2 **h2Data, int dataLen, int NUMPERTURBED);


#endif /* HEADERS_TOFILE_H_ */
