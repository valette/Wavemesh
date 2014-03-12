#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "io.h"

int main( int argc, char *argv[] )
{
	int quantization = 8;
	char *outputFile = 0;

	if(argc < 1)
	{
		std::cout << "Usage :" << std::endl;
		std::cout << "decodeCoeffs inputfile" << std::endl;	
		return (0);
	}

	char* verificationFile = 0;

	// Parse optionnal arguments
	int ArgumentsIndex = 2;
	while (ArgumentsIndex < argc) {
		if (strcmp(argv[ArgumentsIndex],"-v")==0) {
			verificationFile = argv[ArgumentsIndex + 1];
		}
		if (strcmp(argv[ArgumentsIndex],"-o")==0) {
			outputFile = argv[ArgumentsIndex + 1];
		}
		ArgumentsIndex+=2;
	}

	coeffs result = decode(argv[1]);
	if (outputFile) {
		writeFile(result, outputFile);
	} else {
		writeFile(result, "decoded.txt");
	}

	if (verificationFile) {
		coeffs reference = readFile(verificationFile);
		double maxError = 0;
		for (int i = 0; i != reference.size(); i++) {
			coefficient tempCoeff = reference[i];
			for (int j = 0; j != tempCoeff.size(); j++) {
				double error = fabs(tempCoeff[j] - result[i][j]);
				if (maxError < error) {
					maxError = error;
				}
			}
		}
		std::cout << "Maximum error : " << maxError << std::endl;
	}
}

