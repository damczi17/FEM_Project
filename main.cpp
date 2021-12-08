#include<iostream>
#include <algorithm>

#include "net.h"
#include "element4.h"
#include "jakobian.h"
#include "harrays.h"
#include "globalstruct.h"


#define H 0.1 //Height of physical net
#define B 0.1 //Width of physical net
#define nH 4 //Height of net
#define nB 4 //Width of net
#define N 3 //Number of integration points
#define conductivity 25 //Conductivity [W/(mC)]
#define alfa 300 //Conductivity [W/(mC)]
#define ambTemp 1200 //ambient temperature [C]
#define c 700 //specific heat [J/(kgC)]
#define ro 7800 //density [kg/m3]
#define step 50 //simulation step time [s]
#define time 500 //simulation time [s]
#define initTemperature 100 //initial temperature [C]


int main() {

	grid net = netGenerate(H, B, nH, nB, initTemperature);
	element4_2D elem = derivate(N);

	jakobianCnt(net, elem);

	harraycnt(net, elem, conductivity);

	wallsCnt(net, elem, alfa, ambTemp);

	HbcCnt(net, elem);

	globalMatices global1;

	cmatrixCnt(net, elem, ro, c);

	agregation(net, global1.Hglobal, global1.Pglobal, global1.Cglobal);

	//for (int i = 0; i < net.nE; ++i) {
	//	std::cout << "Element: " << i + 1 << std::endl;
	//	for (int j = 0; j < 4; ++j) {
	//		for (int k = 0; k < 4; ++k) {
	//			std::cout << net.elements[i].Hbc[j][k] << " ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	std::cout << std::endl;
	//}
	
	std::cout << net.elements[0].jak.det[0];
	global1.showMatrices();

	globalH_C(net, global1.Hglobal, global1.Cglobal, step);

	//std::cout << "\n\nGlobal H + C/dT matrix\n";
	//for (int i = 0; i < net.nN; ++i) {
	//	for (int j = 0; j < net.nN; ++j) {
	//		std::cout << global1.Hglobal[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}


	std::vector<double> tempSolution;

	for (int i = step; i <= time ; i+=step) {
		std::cout << "\n--------------------------------------------\n";
		
		tempSolution = PvectorCnt(net, global1, step);
		tempSolution = equationSolve(global1.Hglobal, tempSolution);

		auto [min, max] = std::minmax_element(begin(tempSolution), end(tempSolution));
		std::cout << "\nTime " << i << " sec temperature value: " << *min << " Max temperature value: " << *max << std::endl;
		updateNodesTemp(net, tempSolution);
	}

	return 0;
}
