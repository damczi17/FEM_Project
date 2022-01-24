#pragma once
#ifndef import_h
#define import_h
#include<vector>
#include <string>
#include <fstream>
#include <sstream>  

struct Simulation {
	int SimulationTime;//simulation time [s]
	int SimulationStepTime; //simulation step time [s]
	int Conductivity; //Conductivity [W/(mC)]
	int Alfa; //Conductivity [W/(mC)]
	int Tot; //ambient temperature [C]
	int InitialTemp; //initial temperature [C]
	int Density; //specific heat [J/(kgC)]
	int SpecificHeat; //density [kg/m3]
}Simulation1;

void importData(std::string filename, grid& net) {
	
	std::ifstream odczyt;
	odczyt.open(filename);

	std::string a;

	if (odczyt.good()) {
		odczyt >> a >> Simulation1.SimulationTime;
		odczyt >> a >> Simulation1.SimulationStepTime;
		odczyt >> a >> Simulation1.Conductivity;
		odczyt >> a >> Simulation1.Alfa;
		odczyt >> a >> Simulation1.Tot;
		odczyt >> a >> Simulation1.InitialTemp;
		odczyt >> a >> Simulation1.Density;
		odczyt >> a >> Simulation1.SpecificHeat;


		odczyt >> a >> a >> net.nN;
		odczyt >> a >> a >> net.nE;
		odczyt >> a >> a;


		net.elements.resize(net.nE);
		net.nodes.resize(net.nN);

		for (int i = 0; i < net.nN; ++i) {
			net.nodes[i].Temp = Simulation1.InitialTemp; //przypisanie temperatury poczatkowej do wezlow
		}

		std::string line, x1, y1;
		int  i = 0;
		double x2, y2;

		while (getline(odczyt, line)) {
			odczyt >> a;
			std::stringstream ss(line);
			getline(ss, x1, ',');
			x2 = stold(x1);
			getline(ss, y1, ',');
			y2 = stold(y1);
			net.nodes[i].x = x2;
			net.nodes[i].y = y2;
			i++;
			if (i == net.nN)
				break;
		}

		i = 0;
		odczyt >> a >> a;

		std::string id1s, id2s, id3s, id4s;
		int id1, id2, id3, id4;
		while (getline(odczyt, line)) {
			odczyt >> a;
			std::stringstream ss(line);
			getline(ss, id1s, ',');
			id1 = stoi(id1s);
			getline(ss, id2s, ',');
			id2 = stoi(id2s);
			getline(ss, id3s, ',');
			id3 = stoi(id3s);
			getline(ss, id4s, ',');
			id4 = stoi(id4s);

			net.elements[i].ID[0] = id1;
			net.elements[i].ID[1] = id2;
			net.elements[i].ID[2] = id3;
			net.elements[i].ID[3] = id4;
			i++;
			if (i == net.nE)
				break;
		}

		getline(odczyt, line);
		getline(odczyt, line);
		std::stringstream ss(line);
		std::string BCs = "";
		int BC;
		while (getline(ss, BCs, ',')) {
			BC = stoi(BCs);
			net.nodes[BC - 1].BC = 1;
		}


	}
	else {
		std::cout << "Can not read the file!\n";
	}
	odczyt.close();


}

#endif