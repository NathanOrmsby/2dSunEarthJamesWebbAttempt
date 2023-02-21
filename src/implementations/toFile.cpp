/*
 * toFile.cpp
 *
 *  Created on: Feb 19, 2023
 *      Author: norms
 */

#include "../headers/toFile.h"

// Writes h2 probe evolution relative to james webb
void h2Evolution(vector2 *jwData, vector2 **h2Data, int dataLen)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\singleProbeData\\";

	// File name
	std::string f = "h2Evolution.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header
	file << "hx0,hy0,hy1,hy2" << std::endl;

	// Add data to the file
	for (int i = 0; i < dataLen; ++i)
	{
		// h2 positions relative to probe
		double relx0 = h2Data[i][0].x - jwData[i].x;
		double rely0 = h2Data[i][0].y - jwData[i].y;
		double relx1 = h2Data[i][1].x - jwData[i].x;
		double rely1 = h2Data[i][1].y - jwData[i].y;


		file << std::to_string(relx0) << "," << std::to_string(rely0) << "," << std::to_string(relx1) << "," << std::to_string(rely1) << std::endl;
	}

	// Close the file
	file.close();
}
// Writes evolution of circle of probes into ellipse
void perturbedEvolution(vector2 *jwData, vector2 **perturbedData, int dataLen, int NUMPERTURBED)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\singleProbeData\\";

	// File name
	std::string f = "perturbedEvolution.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header
	for (int i = 0; i < NUMPERTURBED - 1; ++i)
	{
		file << "x" + std::to_string(i) + "," + "y" + std::to_string(i) << ",";
	}
	// Last probe: end on newline
	file << "x" + std::to_string(NUMPERTURBED - 1) + "," + "y" + std::to_string(NUMPERTURBED - 1) + "\n";

	// Add data to the file
	for (int i = 0; i < dataLen; ++i)
	{
		// Write to the file
		for (int j = 0; j < NUMPERTURBED - 1; ++j)
		{
			// Stringstream
			std::stringstream x;
			std::stringstream y;

			// Set precision
			x.precision(15);
			y.precision(15);

			// Feed into stringstream
			x << perturbedData[i][j].x - jwData[i].x;
			y << perturbedData[i][j].y - jwData[i].y;

			// Form the string
			std::string out = x.str() + "," + y.str() + ",";

			// Write to file
			file << out;
		}

		// Last perturbed probe: End with newline
		std::stringstream x;
		std::stringstream y;

		x.precision(15);
		y.precision(15);

		x << perturbedData[i][NUMPERTURBED - 1].x - jwData[i].x;
		y << perturbedData[i][NUMPERTURBED - 1].y - jwData[i].y;

		std::string out = x.str() + "," + y.str() + "\n";

		file << out;
	}

	// Close the file
	file.close();
}
// Writes position data of perturbed probes relative to the james webb to csv. This function only works if there are four perturbed probes
void perturbedSeparation(vector2 *jwData, vector2 **perturbedData, int dataLen, int NUMPERTURBED)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\singleProbeData\\";

	// File name
	std::string f = "perturbedSeparation.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File headers
	file << "rightx,righty,topx,topy,bottomx,bottomy,leftx,lefty" << std::endl;

	// Add data to the file
	for (int i = 0; i < dataLen; ++i)
	{
		// Perturbed probes
		// Right
		double rightXRel = perturbedData[i][0].x - jwData[i].x;
		double rightYRel = perturbedData[i][0].y - jwData[i].y;

		// Top
		double topXRel = perturbedData[i][1].x - jwData[i].x;
		double topYRel = perturbedData[i][1].y - jwData[i].y;

		// Left
		double leftXRel = perturbedData[i][2].x - jwData[i].x;
		double leftYRel = perturbedData[i][2].y - jwData[i].y;

		// Bottom
		double bottomXRel = perturbedData[i][3].x - jwData[i].x;
		double bottomYRel = perturbedData[i][3].y - jwData[i].y;

		// Stringstreams
		std::stringstream rightX;
		std::stringstream rightY;
		std::stringstream topX;
		std::stringstream topY;
		std::stringstream bottomX;
		std::stringstream bottomY;
		std::stringstream leftX;
		std::stringstream leftY;

		// Set precision to 15 digits
		rightX.precision(15);
		rightY.precision(15);
		topX.precision(15);
		topY.precision(15);
		bottomX.precision(15);
		bottomY.precision(15);
		leftX.precision(15);
		leftY.precision(15);

		// Feed values into the string streams
		rightX << rightXRel;
		rightY << rightYRel;
		topX << topXRel;
		topY << topYRel;
		bottomX << bottomXRel;
		bottomY << bottomYRel;
		leftX << leftXRel;
		leftY << leftYRel;

		// Form the string
		std::string out = rightX.str() + "," + rightY.str() + "," + topX.str() + "," + topY.str() + "," + bottomX.str() + "," + bottomY.str() + "," + leftX.str() + "," + leftY.str() + "\n";
		file << out;
	}

	// Close the file
	file.close();
}
// Writes position data to csv file for earth and james webb
void earthJWdataToFile(vector2 *earthData, vector2 *jwData, vector2 **perturbedData, vector2 **h2Data, int dataLen, int NUMPERTURBED)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\singleProbeData\\";

	// File name
	std::string f = "earthAndprobe.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File headers
	file << "earthx,earthy,jwx,jwy" << std::endl;

	// TODO: ADD THE DISTANCE DATA TO THE FILE
	for (int i = 0; i < dataLen; ++i)
	{
		// File data:
		// Earth
		file << std::to_string(earthData[i].x) << "," << std::to_string(earthData[i].y) << ",";
		// James Webb
		file << std::to_string(jwData[i].x) << "," << std::to_string(jwData[i].y) << std::endl;
	}

	// Close the file
	file.close();
}


