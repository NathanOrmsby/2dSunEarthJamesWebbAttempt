/*
 * circleGenerator.cpp
 *
 *  Created on: Feb 8, 2023
 *      Author: norms
 */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "../headers/circleGenerator.h"

void writeCircleToFile(int numPoints, double radius)
{
	double **data = (double **)malloc(numPoints * sizeof(double *));
	for (int i = 0; i < numPoints; ++i) { data[i] = (double *)malloc(2 * sizeof(double)); }
	uniformCircle(numPoints, data, radius);

	std::string fname = "unitCircle.csv";
	toFile(fname, data, numPoints);

	// Free data
	for (int i = 0; i < numPoints; ++i) { free(data[i]); } free(data);
}

void uniformCircle(int numPoints, double **arr, double r)
{
	double d = (2.0 * M_PI) / (double)numPoints;
	for (int i = 0; i < numPoints; ++i)
	{
		double theta = i * d;
		arr[i][0] = r * cos(theta); arr[i][2] = r * sin(theta);
	}
}

void toFile(std::string fileName, double **data, int dataLen)
{
	std::ofstream file;
	file.open(fileName);

//	file << "x" << "," << "y" << "," << "z" << std::endl;
	for (int i = 0; i < dataLen; ++i)
	{
		file << std::to_string(data[i][0]) << "," << std::to_string(data[i][1]) << std::endl;
	}
	file.close();
}


