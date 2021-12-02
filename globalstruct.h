#pragma once
#ifndef globalstruct_h
#define globalstruct_h
#include<vector>

std::vector<double> equationSolve(std::vector<std::vector<double>> a, std::vector<double> b) {

	int num = a.size();
	std::vector<double> n, x1, x2;
	std::vector<std::vector<double>> M;
	M.resize(num);
	n.resize(num);
	x1.resize(num);
	x2.resize(num);

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

	for (int k = 0; k < 100000; k++) {
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
#endif
