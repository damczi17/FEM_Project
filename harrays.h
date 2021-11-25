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
	harrays() {
		tabX.resize(4);
		tabY.resize(4);
		for (int i = 0; i < 4; ++i) {
			tabX[i].resize(4);
			tabY[i].resize(4);
		}
	}
};

void harraycnt(grid& net, element4_2D elem, double conductivity) {
	std::vector<harrays> Harray1;

	Harray1.resize(net.nE);
	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				Harray1[i].tabX[j][k] = net.elements[i].jak.oPC[0][0][0] * elem.dKsi[j][k] + net.elements[i].jak.oPC[0][0][1] * elem.dEta[j][k];
				Harray1[i].tabY[j][k] = net.elements[i].jak.oPC[0][1][0] * elem.dKsi[j][k] + net.elements[i].jak.oPC[0][1][1] * elem.dEta[j][k];
			}
		}
	}

		//std::cout << "\nMacierz X\n";

		//for (int i = 0; i < net.nE; ++i) {
		//	for (int j = 0; j < 4; ++j) {
		//		for (int k = 0; k < 4; ++k) {
		//			std::cout << Harray1[i].tabX[j][k] << " ";
		//		}
		//		std::cout << std::endl;
		//	}
		//	std::cout << std::endl;
		//}

		//std::cout << "\nMacierz Y\n";

		//for (int i = 0; i < net.nE; ++i) {
		//	for (int j = 0; j < 4; ++j) {
		//		for (int k = 0; k < 4; ++k) {
		//			std::cout << Harray1[i].tabY[j][k] << " ";
		//		}
		//		std::cout << std::endl;
		//	}
		//	std::cout << std::endl;
		//}

	std::vector<std::vector<vec2D>> tabH, tabH2, sum;
	vec2D result;
	tabH.resize(net.nE);
	tabH2.resize(net.nE);
	sum.resize(net.nE);

	for (int i = 0; i < net.nE; ++i) {
		tabH[i].resize(4);
		tabH2[i].resize(4);
		sum[i].resize(4);
		for (int j = 0; j < 4; ++j) {
			tabH[i][j].resize(4);
			tabH2[i][j].resize(4);
			sum[i][j].resize(4);
			for (int k = 0; k < 4; ++k) {
				tabH[i][j][k].resize(4);
				tabH2[i][j][k].resize(4);
				sum[i][j][k].resize(4);
			}
		}
	}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					tabH[i][j][k][l] = Harray1[i].tabX[j][k] * Harray1[i].tabX[j][l];
					tabH2[i][j][k][l] = Harray1[i].tabY[j][k] * Harray1[i].tabY[j][l];
				}
			}
		}
	}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					sum[i][j][k][l] = conductivity * (tabH[i][j][k][l] + tabH2[i][j][k][l]) * net.elements[i].jak.det[k]; //Wyliczanie macierzy H dla PC
				}
			}
		}
	}

	//for (int i = 0; i < net.nE; ++i) {
	//	for (int j = 0; j < 4; ++j) {
	//		std::cout << "\nH matrix PC " << j + 1 << "\n";	
	//		for (int k = 0; k < 4; ++k) {
	//			for (int l = 0; l < 4; ++l) {
	//				std::cout << sum[i][j][k][l] << " ";
	//			}
	//			std::cout << std::endl;
	//		}
	//	}
	//}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					net.elements[i].Hmatrix[k][l] += sum[i][j][k][l];//Wyliczanie macierzy H dla elementow skoczonych
				}
			}
		}
	}
}

void HbcCnt(grid& net, element4_2D elem) {
	for (int i = 0; i < net.nE; ++i) {
		if (net.elements[i].cords[0].BC == 0 && net.elements[i].cords[2].BC == 0) {
			continue;
		}
		else {
			if (net.elements[i].cords[0].BC == 1 && net.elements[i].cords[1].BC == 1) {
				if (net.elements[i].cords[2].BC == 1) {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[0][j] + elem.Pvector[1][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[0][j][k] + elem.surface[1][j][k];
						}
					}
				}
				else if (net.elements[i].cords[3].BC == 1) {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[0][j] + elem.Pvector[3][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[0][j][k] + elem.surface[3][j][k];
						}
					}
				}
				else {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[0][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[0][j][k];
						}
					}
				}
			}
			else if (net.elements[i].cords[2].BC == 1 && net.elements[i].cords[3].BC == 1) {
				if (net.elements[i].cords[1].BC == 1) {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[1][j] + elem.Pvector[2][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[1][j][k] + elem.surface[2][j][k];
						}
					}
				}
				else if (net.elements[i].cords[0].BC == 1) {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[2][j] + elem.Pvector[3][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[2][j][k] + elem.surface[3][j][k];
						}
					}
				}
				else {
					for (int j = 0; j < 4; ++j) {
						net.elements[i].P[j] = elem.Pvector[2][j];
						for (int k = 0; k < 4; ++k) {
							net.elements[i].Hbc[j][k] = elem.surface[2][j][k];
						}
					}
				}
			}
			else if (net.elements[i].cords[1].BC == 1 && net.elements[i].cords[2].BC == 1) {
				for (int j = 0; j < 4; ++j) {
					net.elements[i].P[j] = elem.Pvector[1][j];
					for (int k = 0; k < 4; ++k) {
						net.elements[i].Hbc[j][k] = elem.surface[1][j][k];
					}
				}
			}
			else if (net.elements[i].cords[0].BC == 1 && net.elements[i].cords[3].BC == 1) {
				for (int j = 0; j < 4; ++j) {
					net.elements[i].P[j] = elem.Pvector[3][j];
					for (int k = 0; k < 4; ++k) {
						net.elements[i].Hbc[j][k] = elem.surface[3][j][k];
					}
				}
			}
		}
	}
}

void globalH(grid net, std::vector<std::vector<double>> &Hglobal, std::vector<double> &Pvector) {
	Hglobal.resize(net.nN);
	Pvector.resize(net.nN);
	for (int i = 0; i < net.nN; ++i) {
		Hglobal[i].resize(net.nN);
	}

	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < 4; ++j) {
			Pvector[net.elements[i].ID[j] - 1] += net.elements[i].P[j];
			for (int k = 0; k < 4; ++k) {
				Hglobal[net.elements[i].ID[j] - 1][net.elements[i].ID[k] - 1] += net.elements[i].Hmatrix[j][k];
			}
		}
	}
}

#endif