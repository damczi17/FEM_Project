#include<iostream>

#include "net.h"
#include "element4.h"
#include "jakobian.h"
#include "harrays.h"


#define H 0.1 //Height of physical net
#define B 0.1 //Width of physical net
#define nH 4 //Height of net
#define nB 4//Width of net
#define N 2 //Number of integration points
#define conductivity 25 //Conductivity [W/(mC)]
#define ambTemp 1200 //ambient temperature [C]

int main() {

	grid net = netGenerate(H, B, nH, nB);
	element4_2D elem = derivate(N);

	jakobianCnt(net, elem);

	harraycnt(net, elem, conductivity);

	wallsCnt(net, elem, conductivity, ambTemp);

	HbcCnt(net, elem);

	//for (int i = 0; i < net.nE; ++i) {
	//	std::cout << "Element " << i + 1 << std::endl;
	//	for (int j = 0; j < 4; ++j) {
	//		for (int k = 0; k < 4; ++k) {
	//			std::cout << net.elements[i].Hbc[j][k] << " ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	std::cout << std::endl;
	//}

	std::vector<std::vector<double>> Hglobal;
	std::vector<double> Pglobal;

	globalH(net, Hglobal, Pglobal);

	std::cout << "\nMacierz Hbc dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << net.elements[0].Hbc[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\nMacierz H dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << net.elements[0].Hmatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\nWektor P dla PC1:\n";
	for (int i = 0; i < 4; ++i) {
		std::cout << net.elements[0].P[i] << "\n";
	}


	std::cout << "\n\nGlobal H matrix\n";
	for (int i = 0; i < net.nN; ++i) {
		for (int j = 0; j < net.nN; ++j) {
			std::cout << Hglobal[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "\n\nGlobal P vector\n";
	for (int i = 0; i < net.nN; ++i) {
		std::cout << Pglobal[i] << "\n ";
	}

	//temperatura musi wyjsc 1200

	return 0;
}
