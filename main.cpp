#include <iostream>
#include <algorithm>

#include "net.h"
#include "element4.h"
#include "jakobian.h"
#include "harrays.h"
#include "globalstruct.h"
#include "import.h"



#define H 0.1 //Height of physical net
#define B 0.1 //Width of physical net
#define nH 5 //Height of net
#define nB 5 //Width of net
#define N 2 //Number of integration points


int main() {

	//ISO C++17 Standard

	//Initial data:
	//Version without imput file:
	
	Simulation1.SimulationTime = 500;//simulation time [s]
	Simulation1.SimulationStepTime = 50; //simulation step time [s]
	Simulation1.Conductivity = 25; //Conductivity [W/(mC)]
	Simulation1.Alfa = 300; //Alfa [W/m^2K)]
	Simulation1.Tot = 1200; //ambient temperature [C]
	Simulation1.InitialTemp = 100; //initial temperature [C]
	Simulation1.Density = 7800;//density [kg/m3]
	Simulation1.SpecificHeat = 700;  //specific heat [J/(kgC)]
	grid net = netGenerate(H, B, nH, nB, Simulation1.InitialTemp);
	
	//end


	//Initial data:
	//Version with input file
	
	//grid net;
	//importData("mynet.txt", net);

	//for (int i = 0; i < net.nE; ++i) {
	//	for (int j = 0; j < 4; ++j) {
	//		net.elements[i].cords[j] = net.nodes[net.elements[i].ID[j] - 1];
	//	}
	//}
	
	//end

	element4_2D elem = derivate(N);

	jakobianCnt(net, elem);

	harraycnt(net, elem, Simulation1.Conductivity);

	wallsCnt(net, elem, Simulation1.Alfa, Simulation1.Tot);

	HbcCnt(net, elem);

	cmatrixCnt(net, elem, Simulation1.Density, Simulation1.SpecificHeat);

	agregation(net, global1.Hglobal, global1.Pglobal, global1.Cglobal, global1.Hbcglobal);
	
	global1.showMatrices();
	//simulation
	globalH_C(net, global1.Hglobal, global1.Cglobal, Simulation1.SimulationStepTime);
	

	global1.showMatrices();

	std::vector<double> tempSolution;

	for (int i = Simulation1.SimulationStepTime; i <= Simulation1.SimulationTime; i+= Simulation1.SimulationStepTime) {
		
		tempSolution = PvectorCnt(net, global1, Simulation1.SimulationStepTime);
		tempSolution = equationSolve(global1.Hglobal, tempSolution);

		for (int i = 0; i < net.nN; ++i) {
			std::cout << tempSolution[i] << std::endl;
		}

		auto [min, max] = std::minmax_element(begin(tempSolution), end(tempSolution));
		std::cout << "\nTime " << i << " sec temperature value: " << *min << " Max temperature value: " << *max << std::endl;
		updateNodesTemp(net, tempSolution);

	}

	return 0;
}
