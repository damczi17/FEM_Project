#pragma once
#ifndef net_h
#define net_h
#include <iostream>
#include <iomanip>
#include <vector>

typedef std::vector<std::vector<double>> vec2D;

struct jakobian {
	std::vector<vec2D> PC, oPC;
	std::vector<double> det;
	jakobian() {
		PC.resize(4);
		oPC.resize(4);
		det.resize(4);
		for (int i = 0; i < 4; ++i) {
			PC[i].resize(2);
			oPC[i].resize(2);
			for (int j = 0; j < 2; ++j) {
				PC[i][j].resize(2);
				oPC[i][j].resize(2);
			}
		}
	}
	void showJakobian() {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 2; ++j) {
				for (int k = 0; k < 2; ++k) {
					std::cout << this->oPC[i][j][k] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "Wyznacznik: " << det[i] << std::endl;
		}
	}
};

struct node {
	double x, y;
	bool BC;
	double initTemp;
	void show() {
		std::cout << "(" << x << "," << y << ") BC: " << BC;
	}
};

struct element {
	int ID[4];
	node cords[4];
	jakobian jak;
	double Hmatrix[4][4];
	double Cmatrix[4][4];
	double Hbc[4][4];
	double P[4];

	void show() {
		std::cout << "( ";
		for (int i = 0; i < 4; ++i)
			std::cout << ID[i] << " ";
		std::cout << ")";
	}
};

struct grid {
	double H, B;  //Net physical dimensions
	int nH, nB; //Nodes net dimensions

	int nN, nE; //nN - number of nodes, nE - number of elements
	double dX, dY; //dX - single step on X axis, dY - single step on Y axis

	std::vector <element> elements;
	std::vector <node> nodes;

	grid(double H, double B, int nH, int nB) {
		this->H = H;
		this->B = B;
		this->nH = nH;
		this->nB = nB;

		this->nN = nH * nB;
		this->nE = (nH - 1) * (nB - 1);
		this->dX = B / (nB - 1);
		this->dY = H / (nH - 1);
	}

	void showInfo() {
		std::cout << "H: " << this->H << "\nB: " << this->B << "\nnH: " << this->nH << "\nnB: " << this->nB << "\nnN: " << this->nN << "\nnE: "
			<< this->nE << "\ndX: " << this->dX << "\ndY: " << this->dY << "\n\n";
	}

	void showNet() {

		std::cout << "Nodes:\n";
		for (int i = this->nH; i >= 1; --i) {
			for (int j = i; j <= this->nN; j += this->nH) {
				nodes[j - 1].show();
				std::cout << "\t\t";
			}
			std::cout << std::endl;
		}

		std::cout << "\n\n";

		std::cout << "Elements:\n";
		for (int i = this->nB; i >= 1; --i) {
			for (int j = i; j <= this->nE; j += this->nB) {
				elements[j - 1].show();
				std::cout << "\t\t";
			}
			std::cout << std::endl;
		}
	}
};



void ID_cnt(int tmp1, int nH, std::vector <element>& elements) {
	int tmp2 = tmp1 + nH;
	elements.push_back({ tmp1, tmp2, tmp2 + 1, tmp1 + 1 });
}

void nodesID(grid &net) {
	double dX = net.dX, dY = net.dY;

	for (int i = 0; i < net.nB; ++i) {
		for (int j = 0; j < net.nH; ++j) {
			if (i == 0 || j == 0 || i + 1 == net.nB || j + 1 == net.nH)
				net.nodes.push_back({ i * dX, j * dY, 1 });
			else
				net.nodes.push_back({ i * dX, j * dY, 0 });
		}
	}

}

void elementsID(grid &net) {
	int tmp = 1, nH = net.nH;
	ID_cnt(tmp, nH, net.elements);
	++tmp;

	for (int i = 1; i < net.nE; ++i) {
		if (tmp % net.nH != 0) {
			ID_cnt(tmp, nH, net.elements);
			++tmp;
		}
		else {
			++tmp;
			ID_cnt(tmp, nH, net.elements);
			++tmp;
		}
	}
}

grid netGenerate(double H, double B, int nH, int nB) {

	grid net = { H, B, nH, nB };

	nodesID(net);

	elementsID(net);

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			net.elements[i].cords[j] = net.nodes[net.elements[i].ID[j] - 1];
		}
	}

	net.showNet();

	return net;
}
#endif
