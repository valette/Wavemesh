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
		std::cout << "encode inputfile" << std::endl;	
		return (0);
	}

	// Parse optionnal arguments
	int ArgumentsIndex = 2;
	while (ArgumentsIndex < argc) {
		if (strcmp(argv[ArgumentsIndex],"-q")==0) {
			quantization = atoi(argv[ArgumentsIndex + 1]);
		}
		if (strcmp(argv[ArgumentsIndex],"-o")==0) {
			outputFile = argv[ArgumentsIndex + 1];
		}
		ArgumentsIndex+=2;
	}

	std::cout << "Quantization : " << quantization << " bits" << std::endl;
	coeffs input = readFile(argv[1]);
	std::cout << input.size() << " coefficients" << std::endl;
	std::cout << input[0].size() << " scalars per coefficient" << std::endl;
	if (outputFile) {
		encode(input, outputFile, quantization);
	} else {
		encode(input, "compressed.svc", quantization);
	}
}

