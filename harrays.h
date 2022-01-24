#pragma once
#ifndef harrays_h
#define harrays_h
#include <iostream>
#include <iomanip>
#include <vector>

#include "net.h"
#include "element4.h"


typedef std::vector<std::vector<double>> vec2D;

struct harrays {
	vec2D tabX, tabY;
	void chng_size(int points) {
		tabX.resize(points);
		tabY.resize(points);
		for (int i = 0; i < points; ++i) {
			tabX[i].resize(4);
			tabY[i].resize(4);
		}
	}
};

//Wyliczanie dN/dx oraz dN/dy
void harraycnt(grid& net, element4_2D elem, double conductivity) {
	
	std::vector<harrays> Harray1;
	Harray1.resize(net.nE);

	for (int i = 0; i < net.nE; ++i) {
		Harray1[i].chng_size(elem.pointsNumber);
		for (int j = 0; j < elem.pointsNumber; ++j) {
			for (int k = 0; k < 4; ++k) {
				Harray1[i].tabX[j][k] = net.elements[i].jak.oPC[j][0][0] * elem.dKsi[j][k] + net.elements[i].jak.oPC[j][0][1] * elem.dEta[j][k];
				Harray1[i].tabY[j][k] = net.elements[i].jak.oPC[j][1][0] * elem.dKsi[j][k] + net.elements[i].jak.oPC[j][1][1] * elem.dEta[j][k];
			}
		}
	}

	std::vector<std::vector<vec2D>> tabH, tabH2;
	tabH.resize(net.nE);
	tabH2.resize(net.nE);

	for (int i = 0; i < net.nE; ++i) {
		tabH[i].resize(elem.pointsNumber);
		tabH2[i].resize(elem.pointsNumber);
		for (int j = 0; j < elem.pointsNumber; ++j) {
			tabH[i][j].resize(4);
			tabH2[i][j].resize(4);
			for (int k = 0; k < 4; ++k) {
				tabH[i][j][k].resize(4);
				tabH2[i][j][k].resize(4);
			}
		}
	}

	//Wyliczanie macierzy H dla punktow calkowania
	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < elem.pointsNumber; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					tabH[i][j][k][l] = Harray1[i].tabX[j][k] * Harray1[i].tabX[j][l];
					tabH2[i][j][k][l] = Harray1[i].tabY[j][k] * Harray1[i].tabY[j][l];
				}
			}
		}
	}


	//Wyliczanie macierzy H dla elementu
	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < elem.pointsNumber; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					net.elements[i].Hmatrix[k][l] += conductivity * elem.ak2[j] * (tabH[i][j][k][l] + tabH2[i][j][k][l]) * net.elements[i].jak.det[k];
				}
			}
		}
	}
}


void HbcMatrixAdd(element &element, vec2D surface, std::vector<double> Pvec) {
	for (int j = 0; j < 4; ++j) {
		element.P[j] += Pvec[j];
		for (int k = 0; k < 4; ++k) {
			element.Hbc[j][k] += surface[j][k];
		}
	}
}
//Wyliczanie macierzy Hbc i wketroa P dla elementu
void HbcCnt(grid& net, element4_2D elem) {
	for (int i = 0; i < net.nE; ++i) {
		if (net.elements[i].cords[0].BC == 1 && net.elements[i].cords[1].BC == 1) {
			HbcMatrixAdd(net.elements[i], elem.surface[0], elem.Pvector[0]);
		}
		if (net.elements[i].cords[1].BC == 1 && net.elements[i].cords[2].BC == 1) {
			HbcMatrixAdd(net.elements[i], elem.surface[1], elem.Pvector[1]);
		}
		if (net.elements[i].cords[2].BC == 1 && net.elements[i].cords[3].BC == 1) {
			HbcMatrixAdd(net.elements[i], elem.surface[2], elem.Pvector[2]);
		}
		if (net.elements[i].cords[3].BC == 1 && net.elements[i].cords[0].BC == 1) {
			HbcMatrixAdd(net.elements[i], elem.surface[3], elem.Pvector[3]);
		}
	}
}


//Agregacja globalnych macierzy
void agregation(grid net, std::vector<std::vector<double>> &Hglobal, std::vector<double> &Pvector, std::vector<std::vector<double>>& Cglobal, std::vector<std::vector<double>>& Hbcglobal) {
	Hglobal.resize(net.nN);
	Hbcglobal.resize(net.nN);
	Cglobal.resize(net.nN);
	Pvector.resize(net.nN);
	for (int i = 0; i < net.nN; ++i) {
		Hbcglobal[i].resize(net.nN);
		Hglobal[i].resize(net.nN);
		Cglobal[i].resize(net.nN);
	}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			Pvector[net.elements[i].ID[j] - 1] += net.elements[i].P[j];
			for (int k = 0; k < 4; ++k) {
				Hglobal[net.elements[i].ID[j] - 1][net.elements[i].ID[k] - 1] += (net.elements[i].Hmatrix[j][k]);
				Hbcglobal[net.elements[i].ID[j] - 1][net.elements[i].ID[k] - 1] += (net.elements[i].Hbc[j][k]);
				Cglobal[net.elements[i].ID[j] - 1][net.elements[i].ID[k] - 1] += net.elements[i].Cmatrix[j][k];
			}
		}
	}
}

#endif