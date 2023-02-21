/*
 * multiProbeToFile.h
 *
 *  Created on: Feb 20, 2023
 *      Author: norms
 */

#ifndef HEADERS_MULTIPROBETOFILE_H_
#define HEADERS_MULTIPROBETOFILE_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>

#include "vector2.h"

// Writes position data to csv file for earth and james webb
void earthAndProbesToFile(vector2 *earthData, vector2 **probesData, int dataLen, int numProbes);

// Writes position data to csv file for probes, perturbations, and h2 probes
void earthAndProbesToFile(vector2 *earthData, vector2 **probesData, vector2 ***perturbedData, vector2 ***h2Data, int dataLen, int numProbes, int NUMPERTURBED);

// Writes lyapunov exponents to file in width x height format
void lyapunovsToFile(double **lyapunovs, int width, int height);

// All probe perturbation evolutions
void multiProbePerturbationEvolution(vector2 **probesData, vector2 ***perturbedData, int dataLen, int numProbes, int NUMPERTURBED);



#endif /* HEADERS_MULTIPROBETOFILE_H_ */
