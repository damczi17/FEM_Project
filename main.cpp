#include<iostream>

#include "net.h"
#include "element4.h"
#include "jakobian.h"
#include "harrays.h"
#include "globalstruct.h"


#define H 0.1 //Height of physical net
#define B 0.1 //Width of physical net
#define nH 5 //Height of net
#define nB 5//Width of net
#define N 2 //Number of integration points
#define conductivity 25 //Conductivity [W/(mC)]
#define ambTemp 1200 //ambient temperature [C]
#define c 700 //specific heat [J/(kgC)]
#define ro 7800 //density [kg/m3]
#define step 50 //simulation step time [s]
#define initTemp 100 //initial temperature [C]


int main() {

	grid net = netGenerate(H, B, nH, nB);
	element4_2D elem = derivate(N);

	jakobianCnt(net, elem);

	harraycnt(net, elem, conductivity);

	wallsCnt(net, elem, conductivity, ambTemp);

	HbcCnt(net, elem);

	std::vector<std::vector<double>> Hglobal, Cglobal;
	std::vector<double> Pglobal;

	cmatrixCnt(net, elem, ro, c);

	globalH(net, Hglobal, Pglobal, Cglobal);

	std::cout << "\nMacierz H dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << net.elements[0].Hmatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\nMacierz Hbc dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << net.elements[0].Hbc[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\nMacierz C dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << net.elements[0].Cmatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	
	for (int j = 0; j < net.nE; ++j) {
		std::cout << "\nWektor P dla PC" << j+1 <<"\n";
		for (int i = 0; i < 4; ++i) {
			std::cout << net.elements[j].P[i] << "\n";
		}
	}

	std::cout << "\n\nGlobal H + Hbc matrix\n";
	for (int i = 0; i < net.nN; ++i) {
		for (int j = 0; j < net.nN; ++j) {
			std::cout << Hglobal[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n\nGlobal C matrix\n";
	for (int i = 0; i < net.nN; ++i) {
		for (int j = 0; j < net.nN; ++j) {
			std::cout << Cglobal[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n\nGlobal P vector\n";
	for (int i = 0; i < net.nN; ++i) {
		std::cout << Pglobal[i] << "\n ";
	}

	std::vector<double> tempSolution = equationSolve(Hglobal, Pglobal);

	std::cout << "\n\nSolutions of temperature:\n";
	for (int i = 0; i < net.nN; ++i) {
		std::cout << "Temperature " << i+1 << " " << tempSolution[i] << "\n ";
	}

	return 0;
}
