#ifndef element4_h
#define element4_h
#include <iostream>
#include <iomanip>
#include <vector>

#include "net.h"
#include "element4.h"

typedef std::vector<std::vector<double>> vec2D;

double detCnt(vec2D res) {
	return res[0][0] * res[1][1] - (res[1][0] * res[0][1]);
}

void jcnt(int i, int j, element4_2D elem, grid& net) {
		for (int k = 0; k < 4; ++k) {
			net.elements[i].jak.PC[j][0][0] += (double)elem.dEta[j][k] * (double)net.elements[i].cords[k].y;
			net.elements[i].jak.PC[j][0][1] += (-1.)*((double)elem.dKsi[j][k] * (double)net.elements[i].cords[k].y);
			net.elements[i].jak.PC[j][1][0] += (-1.) * ((double)elem.dEta[j][k] * (double)net.elements[i].cords[k].x);
			net.elements[i].jak.PC[j][1][1] += (double)elem.dKsi[j][k] * (double)net.elements[i].cords[k].x;
		}
}

void oJakobian(grid &net, int element, int point) {
	for (int b = 0; b < 2; ++b) {
		for (int c = 0; c < 2; ++c)
			net.elements[element].jak.oPC[point][b][c] = (1. / (double)net.elements[element].jak.det[point]) * net.elements[element].jak.PC[point][b][c];
	}

}

void jakobianCnt(grid &net, element4_2D &elem) {

	for (int i = 0; i < net.nE; ++i) {
		net.elements[i].jak.chng_size(elem.pointsNumber);
		for (int j = 0; j < elem.pointsNumber; ++j) {
			jcnt(i, j, elem, net);
		}
	}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < elem.pointsNumber; ++j) {
			net.elements[i].jak.det[j] = detCnt(net.elements[i].jak.PC[j]);
		}
	}

	std::cout << std::endl;
	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < elem.pointsNumber; ++j) {
			oJakobian(net, i, j);
		}
	}

	//std::cout << "jakobian arays:\n";
	//for (int i = 0; i < net.nE; ++i) {
	//	std::cout << "element " << i + 1 << std::endl << "----------------------------------------------\n";
	//	net.elements[i].jak.showJakobian(elem.pointsNumber);
	//}
}

#endif
