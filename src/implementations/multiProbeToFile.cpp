/*
 * multiProbeToFile.cpp
 *
 *  Created on: Feb 20, 2023
 *      Author: norms
 */


#include "../headers/multiProbeToFile.h"

// All probe perturbation evolutions with probe as well
void multiProbePerturbationEvolution(vector2 **probesData, vector2 ***perturbedData, int dataLen, int numProbes, int NUMPERTURBED)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\multiProbeData\\";

	// File name
	std::string f = "perturbedEvolution.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header

	// x and y for each perturbedij combo
	for (int i = 0; i < numProbes - 1; ++i)
	{
		// For all perturbed probes
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			file << "x" + std::to_string(i) + std::to_string(j) + "," + "y" + std::to_string(i) + std::to_string(j) << ",";
		}
	}

	for (int j = 0; j < NUMPERTURBED - 1; ++j)
	{
		file << "x" + std::to_string(numProbes - 1) + std::to_string(j) + "," + "y" + std::to_string(numProbes - 1) + std::to_string(j) << ",";
	}
	file << "x" + std::to_string(numProbes - 1) + std::to_string(NUMPERTURBED - 1) + "," + "y" + std::to_string(numProbes - 1) + std::to_string(NUMPERTURBED - 1) << "\n";


	// Add data to the file
	for (int i = 0; i < dataLen; ++i)
	{
		// Write to the file
		for (int j = 0; j < numProbes - 1; ++j)
		{
			// Perturbed
			for (int k = 0; k < NUMPERTURBED; ++k)
			{
				file << std::to_string(perturbedData[i][j][k].x - probesData[i][j].x) << "," << std::to_string(perturbedData[i][j][k].y - probesData[i][j].y) << ",";
			}
		}
		for (int j = 0; j < NUMPERTURBED - 1; ++j)
		{
			file << std::to_string(perturbedData[i][numProbes - 1][j].x - probesData[i][numProbes - 1].x) << "," + std::to_string(perturbedData[i][numProbes - 1][j].y - probesData[i][numProbes - 1].y) << ",";
		}
		file << std::to_string(perturbedData[i][numProbes - 1][NUMPERTURBED - 1].x - probesData[i][numProbes - 1].x) << "," << std::to_string(perturbedData[i][numProbes - 1][NUMPERTURBED - 1].y - probesData[i][numProbes - 1].y) << std::endl;
	}

	// Close the file
	file.close();
}
// Writes lyapunov exponents to file in width x height format
void lyapunovsToFile(double **lyapunovs, int width, int height)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\multiProbeData\\";

	// File name
	std::string f = "lyapunovs.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	int numProbes = width * height;

	file << "lmax,lmin" << std::endl;

	// Add the data to file
	for (int i = 0; i < numProbes; ++i)
	{
		// String streams
		std::stringstream lmax;
		std::stringstream lmin;

		// Custom precision
		lmax.precision(20);
		lmin.precision(20);

		// Scale it to 1 / days from 1 / seconds
		double lmaxScaled = lyapunovs[i][0] * 86400.0;
		double lminScaled = lyapunovs[i][1] * 86400.0;

		// Add to string streams
		lmax << lmaxScaled;
		lmin << lminScaled;

		file << lmax.str() << "," << lmin.str() << std::endl;
	}

	// Close the file
	file.close();
}
// Writes position data to csv file for probes, perturbations, and h2 probes
void earthAndProbesToFile(vector2 *earthData, vector2 **probesData, vector2 ***perturbedData, vector2 ***h2Data, int dataLen, int numProbes, int NUMPERTURBED)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\multiProbeData\\";

	// File name
	std::string f = "probesPerturbedAndh2.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header

	for (int i = 0; i < numProbes - 1; ++i)
	{
		// Probes
		file << "probeX" + std::to_string(i) + "," + "probeY" + std::to_string(i) + ",";

		// Perturbed Probes
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			file << "perturbedX" + std::to_string(i) + std::to_string(j) + "," + "perturbedY" + std::to_string(i) + std::to_string(j) + ",";
		}

		// h2 Probes
		for (int j = 0; j < 2; ++j)
		{
			file << "h2X" + std::to_string(i) + std::to_string(j) + "," + "h2Y" + std::to_string(i) + std::to_string(j) + ",";
		}
	}

	// Last one to end the file
	// Probes
	file << "probeX" + std::to_string(numProbes - 1) + "," + "probeY" + std::to_string(numProbes - 1) + ",";

	// Perturbed Probes
	for (int j = 0; j < NUMPERTURBED; ++j)
	{
		file << "perturbedX" + std::to_string(numProbes - 1) + std::to_string(j) + "," + "perturbedY" + std::to_string(numProbes - 1) + std::to_string(j) + ",";
	}

	// h2 Probes
	file << "h2X" + std::to_string(numProbes - 1) + "0" + "," + "h2Y" + std::to_string(numProbes - 1) + "0" + ",";
	file << "h2X" + std::to_string(numProbes - 1) + "1" + "," + "h2Y" + std::to_string(numProbes - 1) + "1" + "\n";



	// Add the data to file
	for (int i = 0; i < dataLen; ++i)
	{
		// File data:

		// Probes:
		for (int j = 0; j < numProbes - 1; ++j)
		{
			// String streams
			std::stringstream probeX;
			std::stringstream probeY;

			// Custom precision
			probeX.precision(15);
			probeY.precision(15);

			// Add to string streams
			probeX << probesData[i][j].x;
			probeY << probesData[i][j].y;

			// Add to the file
			file << probeX.str() << "," << probeY.str() << ",";

			// Perturbed Probes:
			for (int k = 0; k < NUMPERTURBED; ++k)
			{
				// String streams
				std::stringstream perturbedX;
				std::stringstream perturbedY;

				// Custom precision
				perturbedX.precision(15);
				perturbedY.precision(15);

				// Add to string streams
				perturbedX << perturbedData[i][j][k].x;
				perturbedY << perturbedData[i][j][k].y;

				// Add to the file
				file << perturbedX.str() << "," << perturbedY.str() << ",";
			}

			// h2 Probes:
			for (int k = 0; k < 2; ++k)
			{
				// String streams
				std::stringstream h2X;
				std::stringstream h2Y;

				// Custom precision
				h2X.precision(15);
				h2Y.precision(15);

				// Add to string streams
				h2X << h2Data[i][j][k].x;
				h2Y << h2Data[i][j][k].y;

				// Add to the file
				file << h2X.str() << "," << h2Y.str() << ",";
			}
		}

		// Last probe
		// String streams
		std::stringstream probeX;
		std::stringstream probeY;

		// Custom precision
		probeX.precision(15);
		probeY.precision(15);

		// Add to string streams
		probeX << probesData[i][numProbes - 1].x;
		probeY << probesData[i][numProbes - 1].y;

		// Add to the file
		file << probeX.str() << "," << probeY.str() << ",";

		// Perturbed probes
		for (int j = 0; j < NUMPERTURBED; ++j)
		{
			// String streams
			std::stringstream perturbedX;
			std::stringstream perturbedY;

			// Custom precision
			perturbedX.precision(15);
			perturbedY.precision(15);

			// Add to string streams
			perturbedX << perturbedData[i][numProbes - 1][j].x;
			perturbedY << perturbedData[i][numProbes - 1][j].y;

			// Add to the file
			file << perturbedX.str() << "," << perturbedY.str() << ",";
		}

		// h2 Probes:

		// h2 Probe 0:
		// String streams
		std::stringstream h2X;
		std::stringstream h2Y;

		// Custom precision
		h2X.precision(15);
		h2Y.precision(15);

		// Add to string streams
		h2X << h2Data[i][numProbes - 1][0].x;
		h2Y << h2Data[i][numProbes - 1][0].y;

		// Add to the file
		file << h2X.str() << "," << h2Y.str() << ",";

		// h2 Probe 0:
		// String streams
		std::stringstream h2XLast;
		std::stringstream h2YLast;

		// Custom precision
		h2XLast.precision(15);
		h2YLast.precision(15);

		// Add to string streams
		h2XLast << h2Data[i][numProbes - 1][1].x;
		h2YLast << h2Data[i][numProbes - 1][1].y;

		// Add to the file
		file << h2XLast.str() << "," << h2YLast.str() << "\n";
	}

	// Close the file
	file.close();
}

// Writes position data to csv file for earth and james webb
void earthAndProbesToFile(vector2 *earthData, vector2 **probesData, int dataLen, int numProbes)
{
	// File path
	std::string fpath = "D:\\2dSunEarthChaos\\multiProbeData\\";

	// File name
	std::string f = "earthAndProbes.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header
	file << "earthx,earthy,";
	for (int i = 0; i < numProbes - 1; ++i)
	{
		file << "probeX" + std::to_string(i) + "," + "probeY" + std::to_string(i) + ",";
	}
	file << "probeX" + std::to_string(numProbes - 1) + "," + "probeY" + std::to_string(numProbes - 1) + "\n";


	// Add the data to file
	for (int i = 0; i < dataLen; ++i)
	{
		// File data:

		// Earth:
		// String streams
		std::stringstream earthX;
		std::stringstream earthY;

		// Custom precision
		earthX.precision(15);
		earthY.precision(15);

		// Add to string streams
		earthX << earthData[i].x;
		earthY << earthData[i].y;

		// Add to the file
		file << earthX.str() << "," << earthY.str() << ",";

		// Probes:
		for (int j = 0; j < numProbes - 1; ++j)
		{
			// String streams
			std::stringstream probeX;
			std::stringstream probeY;

			// Custom precision
			probeX.precision(15);
			probeY.precision(15);

			// Add to string streams
			probeX << probesData[i][j].x;
			probeY << probesData[i][j].y;

			// Add to the file
			file << probeX.str() << "," << probeY.str() << ",";
		}

		// Last probe to end the line
		// String streams
		std::stringstream probeX;
		std::stringstream probeY;

		// Custom precision
		probeX.precision(15);
		probeY.precision(15);

		// Add to string streams
		probeX << probesData[i][numProbes - 1].x;
		probeY << probesData[i][numProbes - 1].y;

		// Add to the file
		file << probeX.str() << "," << probeY.str() << "\n";
	}

	// Close the file
	file.close();
}

