#include <vector>
#include "vtkArithmeticCoderBase.h"
#include "QuasiStaticModel.h"

typedef std::vector<double> coefficient;
typedef std::vector<coefficient> coeffs;

void encode(coeffs coefficients, const char *outputFile, int quantization) {
	double max = 0;
	for (int i = 0; i != coefficients.size(); i++) {
		coefficient tempCoeff = coefficients[i];
		for (int j = 0; j != tempCoeff.size(); j++) {
			double temp = fabs(tempCoeff[j]);
			if (max < temp) {
				max = temp;
			}
		}
	}

	std::cout<< "maximum absolute scalar value : " << max << std::endl;
	vtkArithmeticCoderBase *Coder = vtkArithmeticCoderBase::New();
	Coder->OpenFile(outputFile, 1);
	Coder->StartCoding();
	Coder->EncodeByte(quantization);
	Coder->EncodeFloat(max);
	Coder->EncodeInt(coefficients.size());
	Coder->EncodeByte(coefficients[0].size());

	double imax = 2 << (quantization - 1) - 1;
	qsmodel model;
	model.initqsmodel(imax + 1, 15, 1000, NULL, 1);
	for (int i = 0; i != coefficients.size(); i++) {
		coefficient tempCoeff = coefficients[i];
		for (int j = 0; j != tempCoeff.size(); j++) {
			double scalar = tempCoeff[j];
			int iScalar = floor (0.5 + imax * fabs(scalar) / max);
		//std::cout<<iScalar<<std::endl;
			Coder->Encode(iScalar, &model);
			if (iScalar != 0) {
				Coder->EncodeBit(scalar<0);
			}
		}
	}


	int size = Coder->StopCoding();
	cout << size << " bytes written" << std::endl;
	Coder->CloseFile();
}

coeffs decode(const char *inputFile) {
	coeffs coefficients;
	vtkArithmeticCoderBase *Coder = vtkArithmeticCoderBase::New();
	Coder->OpenFile(inputFile, 0);
	Coder->StartDecoding();
	int quantization = Coder->DecodeByte();
	double max = Coder->DecodeFloat();
	std::cout << "Quantization : " << quantization << " bits/scalar" << std::endl;
	std::cout << "Max scalar value : " << max << std::endl;

	int numCoeffs = Coder->DecodeInt();
	int coeffSize = Coder->DecodeByte();
	std::cout << numCoeffs << " vectors" << std::endl;
	std::cout << coeffSize << " scalars per vector" << std::endl;

	double imax = 2 << (quantization - 1) - 1;

	qsmodel model;
	model.initqsmodel(imax + 1, 15, 1000, NULL, 0);
	for (int i = 0; i != numCoeffs; i++) {
		coefficient tempCoeff;
		for (int j = 0; j != coeffSize; j++) {
			double iScalar = Coder->Decode(&model);
			double scalar = iScalar * max / (double) imax;
			if (iScalar != 0) {
				if (Coder->DecodeBit()) {
					scalar = -scalar;
				}
			}
//			std::cout << scalar << std::endl;
			tempCoeff.push_back(scalar);
		}
		coefficients.push_back(tempCoeff);
	}
	Coder->StopDecoding();
	return coefficients;
}

void writeFile (coeffs input, const char *file) {
	fstream	output;
	output.open (file, ofstream::out | ofstream::trunc);
	for (int i = 0; i != input.size(); i++) {
		coefficient tempCoeff = input[i];
		for (int j = 0; j != tempCoeff.size(); j++) {
			output << tempCoeff[j];
			if (j != tempCoeff.size() -1) {
				output << " ";
			} else {
				if (i != input.size() -1) {
					output << std::endl;
				}
			}
		}
	}
}


coeffs readFile (const char *file, int vectorSize) {
	coeffs result;
	std::string line;

	std::ifstream pFile (file);

	if (pFile.is_open()) {
		while(!pFile.eof()) {
			coefficient tempCoeff;
			getline(pFile, line);
			std::stringstream ss(line);
			double coeff;
			while(ss >> coeff) {
				tempCoeff.push_back(coeff);
			}
			if (tempCoeff.size()) {
				if (tempCoeff.size() != vectorSize) {
					std::cout << "Error while reading file!" << std::endl;
					std::cout << "only " << tempCoeff.size() << "coefficients where read" << std::endl;
					return result;
				}
				result.push_back(tempCoeff);
			}
//			if (result.size()>10) break;
		} 
		pFile.close();
		return result;
	}
	else {
		std::cout << "Unable to open file";
	}
}
