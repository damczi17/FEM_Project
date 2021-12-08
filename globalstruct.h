#pragma once
#ifndef globalstruct_h
#define globalstruct_h
#include<vector>

struct globalMatices {
	std::vector<std::vector<double>> Hglobal, Cglobal;
	std::vector<double> Pglobal;

	void showMatrices() {
		std::cout << "\n\nGlobal H + Hbc matrix\n";
		for (int i = 0; i < Hglobal.size(); ++i) {
			for (int j = 0; j < Hglobal.size(); ++j) {
				std::cout << Hglobal[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\n\nGlobal C matrix\n";
		for (int i = 0; i < Cglobal.size(); ++i) {
			for (int j = 0; j < Cglobal.size(); ++j) {
				std::cout << Cglobal[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\n\nGlobal P vector\n";
		for (int i = 0; i < Pglobal.size(); ++i) {
			std::cout << Pglobal[i] << "\n ";
		}
		
	}
};

std::vector<double> equationSolve(std::vector<std::vector<double>> a, std::vector<double> b) {

	int num = a.size();
	std::vector<double> n(num), x1(num), x2(num);
	std::vector<std::vector<double>> M(num);

	for (int i = 0; i < num; i++)
		M[i].resize(num);

	for (int i = 0; i < num; i++)
		n[i] = 1 / a[i][i];

	// Calculate M = -D^-1 (L + U)
	for (int i = 0; i < num; i++)
		for (int j = 0; j < num; j++)
			if (i == j)
				M[i][j] = 0;
			else
				M[i][j] = -(a[i][j] * n[i]);

	for (int k = 0; k < 100; k++) {
		for (int i = 0; i < num; i++) {
			x2[i] = n[i] * b[i];
			for (int j = 0; j < num; j++)
				x2[i] += M[i][j] * x1[j];
		}
		for (int i = 0; i < num; i++)
			x1[i] = x2[i];
	}
	return x1;
}

void updateNodesTemp(grid &net, std::vector<double> Tempresults) {
	for (int i = 0; i < net.nN; ++i) {
		net.nodes[i].Temp = Tempresults[i];
	}
}

void globalH_C(grid net, std::vector<std::vector<double>>& Hglobal, std::vector<std::vector<double>> Cglobal, double step) {
	for (int i = 0; i < net.nN; ++i) {
		for (int j = 0; j < net.nN; ++j) {
			Hglobal[i][j] += (Cglobal[i][j] / step);
		}
	}
}

std::vector<double> PvectorCnt(grid net, globalMatices matrices, double step) {
	std::vector<double> tmp;
	tmp.resize(matrices.Pglobal.size());

	for (int i = 0; i < net.nN; ++i) {
		for (int j = 0; j < net.nN; ++j) {
			tmp[i] += (matrices.Cglobal[i][j]/step) * net.nodes[j].Temp;
		}
		tmp[i] += matrices.Pglobal[i];
	}

	return tmp;
}

#endif
