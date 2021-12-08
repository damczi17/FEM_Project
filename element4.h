#ifndef jakobian_h
#define jakobian_h
#include<iostream>
#include<math.h>
#include<vector>

typedef std::vector<std::vector<double>> vec2D;


struct pom {
	vec2D PC1, PC2, PC3;
	std::vector<double> Pvec;
	pom() {
		PC1.resize(4);
		PC2.resize(4);
		PC3.resize(4);
		Pvec.resize(4);
		for (int i = 0; i < 4; ++i) {
			PC1[i].resize(4);
			PC2[i].resize(4);
			PC3[i].resize(4);
		}
	}

	void showPC() {
		std::cout << "\nPC1:\n-----------------------------\n";
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << this->PC1[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\nPC2:\n-----------------------------\n";
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << this->PC2[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}
};

struct element4_2D {
	int N, pointsNumber;
	std::vector<double> wKsi, wEta, Hcords, ak, ak2;
	vec2D dKsi, dEta, Pvector;
	std::vector<vec2D> walls, surface;


	element4_2D(int N) {
		this->N = N;
		this->pointsNumber = N * N;

		ak.resize(N);
		Hcords.resize(N);
		dKsi.resize(pointsNumber);
		dEta.resize(pointsNumber);
		Pvector.resize(4);

		for (int i = 0; i < pointsNumber; ++i) {
			dKsi[i].resize(4);
			dEta[i].resize(4);
		}

		walls.resize(4);
		surface.resize(4);

		for (int i = 0; i < 4; ++i) {
			walls[i].resize(4);
			surface[i].resize(4);
			Pvector[i].resize(4);
			for (int j = 0; j < 4; ++j) {
				walls[i][j].resize(4);
				surface[i][j].resize(4);
			}
		}

		if (N == 2) {
			double val = 1. / sqrt(3);
			this->wKsi = { -1. * val , val, val, -1. * val };//Valuse on local X axis for 2 integral points
			this->wEta = { -1. * val , -1. * val, val, val };//Valuse on local Y axis for 3 integral points
			this->Hcords = { 1, val };
			this->ak = { 1,1 };
			this->ak2 = { 1,1,1,1 };
		}
		else if (N == 3) {
			ak.resize(pointsNumber);
			double val = sqrt(3. / 5.);
			this->wKsi = { -1. * val , 0, val, -1. * val, 0, val, -1. * val, 0, val };
			this->wEta = { -1. * val , -1. * val, -1. * val, 0, 0, 0, val, val, val };
			this->Hcords = { -1. * val, 0, val };
			this->ak = {5./9. , 8. / 9., 5. / 9. };
			this->ak2 = { 5. / 9. * 5. / 9., 8. / 9. * 5. / 9. , 5. / 9. * 5. / 9, 
				8. / 9. * 5. / 9., 8. / 9. * 8. / 9., 8. / 9. * 5./9.,5. / 9. * 5. / 9., 8. / 9. * 5. / 9. , 5. / 9. * 5. / 9};
		}
	}

	void showDerivates() {
		std::cout << "Funkcje ksztaltu dN/dKsi:\n";
		for (int i = 0; i < pointsNumber; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dKsi[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\nFunkcje ksztaltu dN/dEta:\n";
		for (int i = 0; i < pointsNumber; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dEta[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	void showWalls() {
		std::cout << "\nWalls:\n";
		for (int i = 0; i < 4; ++i) {
			std::cout << "\nWall: " << i + 1 << "\n---------------------------------\n";
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					std::cout << this->surface[i][j][k] << " ";
				}
				std::cout << std::endl;
			}
		}
	}
};

double N1(double x, double y) {
	return (0.25) * (1. - x) * (1. - y);
}
double N2(double x, double y) {
	return (0.25) * (1. + x) * (1. - y);
}
double N3(double x, double y) {
	return (0.25) * (1. + x) * (1. + y);
}
double N4(double x, double y) {
	return (0.25) * (1. - x) * (1. + y);
}


double dN1(double x) {
	return (-1) * (0.25) * (1. - x);
}
double dN2(double x) {
	return 0.25 * (1. - x);
}
double dN3(double x) {
	return 0.25 * (1. + x);
}
double dN4(double x) {
	return (-1) * (0.25) * (1. + x);
}

double detjCnt(grid net, int i) {
	double a = net.elements[0].cords[i+1].x - net.elements[0].cords[i].x;
	double b = net.elements[0].cords[i+1].y - net.elements[0].cords[i].y;
	double res = sqrt((a * a) + (b * b));
	return res / 2;
}

void wallsCnt(grid net, element4_2D& elem, double alfa, double ambTemp) {
	double (*pFunc[4])(double, double);//Funkcje ksztaltu
	pFunc[0] = N1;
	pFunc[1] = N2;
	pFunc[2] = N3;
	pFunc[3] = N4;

	std::vector<double> cords;
	if (elem.N == 2) {

		cords.resize(2);

		for (int i = 0; i < 4; ++i) {
			if (i == 0) {
				cords[0] = (-1. * elem.Hcords[1]);
				cords[1] = (-1. * elem.Hcords[0]);
			}
			else if (i == 1) {
				cords[0] = (elem.Hcords[0]);
				cords[1] = (-1. * elem.Hcords[1]);
			}
			else if (i == 2) {
				cords[0] = (elem.Hcords[1]);
				cords[1] = (elem.Hcords[0]);
			}
			else if (i == 3) {
				cords[0] = (-1. * elem.Hcords[0]);
				cords[1] = (elem.Hcords[1]);
			}

			for (int j = 0; j < elem.N; ++j) {
				if (i == 0 && j == 1)
					cords[0] = elem.Hcords[1];
				else if (i == 1 && j == 1)
					cords[1] = elem.Hcords[1];
				else if (i == 2 && j == 1)
					cords[0] = -1. * elem.Hcords[1];
				else if (i == 3 && j == 1)
					cords[1] = -1. * elem.Hcords[1];
				for (int k = 0; k < 4; ++k) {
					elem.walls[i][j][k] = pFunc[k](cords[0], cords[1]);
				}
			}
		}

		std::vector<pom> pom1;
		pom1.resize(4);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				elem.Pvector[i][j] = alfa * ambTemp * (elem.walls[i][0][j] + elem.walls[i][1][j]) * detjCnt(net, i);//Obliczanie wektora P
				std::cout << detjCnt(net, i);
				for (int k = 0; k < 4; ++k) {
					pom1[i].PC1[j][k] = alfa * (elem.walls[i][0][j] * elem.walls[i][0][k]);
					pom1[i].PC2[j][k] = alfa * (elem.walls[i][1][j] * elem.walls[i][1][k]);
				}
			}
		}

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					elem.surface[i][j][k] = detjCnt(net, i) * (pom1[i].PC1[j][k] + pom1[i].PC2[j][k]);//obliczanie wartosci dla macierzy Hbc 
				}
			}
		}
	}

	else if (elem.N == 3) {
		cords.resize(2);

		for (int i = 0; i < 4; ++i) {
			if (i == 0) {
				cords[0] = (elem.Hcords[0]);
				cords[1] = (-1);
			}
			else if (i == 1) {
				cords[0] = (1);
				cords[1] = (elem.Hcords[0]);
			}
			else if (i == 2) {
				cords[0] = (elem.Hcords[2]);
				cords[1] = (1);	
			}
			else if (i == 3) {
				cords[0] = (-1);
				cords[1] = (elem.Hcords[2]);
			}

			for (int j = 0; j < elem.N; ++j) {
				if (i == 0 && j == 1)
					cords[0] = elem.Hcords[1];
				else if (i == 0 && j == 2)
					cords[0] = elem.Hcords[2];
				else if (i == 1 && j == 1)
					cords[1] = elem.Hcords[1];
				else if (i == 1 && j == 2)
					cords[1] = elem.Hcords[2];
				else if (i == 2 && j == 1)
					cords[0] = elem.Hcords[1];
				else if (i == 2 && j == 2)
					cords[0] = elem.Hcords[0];
				else if (i == 3 && j == 1)
					cords[1] = elem.Hcords[1];
				else if (i == 3 && j == 2)
					cords[1] = elem.Hcords[0];
				for (int k = 0; k < 4; ++k) {
					elem.walls[i][j][k] = pFunc[k](cords[0], cords[1]);
				}
			}
		}

		std::vector<pom> pom1;
		pom1.resize(4);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				elem.Pvector[i][j] = alfa * ambTemp * ((elem.ak[0] * elem.walls[i][0][j]) + (elem.ak[1] * elem.walls[i][1][j]) + (elem.ak[2] * elem.walls[i][2][j])) * detjCnt(net, i);
				for (int k = 0; k < 4; ++k) {
					pom1[i].PC1[j][k] = alfa * elem.ak[0] * (elem.walls[i][0][j] * elem.walls[i][0][k]);
					pom1[i].PC2[j][k] = alfa * elem.ak[1] * (elem.walls[i][1][j] * elem.walls[i][1][k]);
					pom1[i].PC3[j][k] = alfa * elem.ak[2] * (elem.walls[i][2][j] * elem.walls[i][2][k]);
				}
			}
		}
		//Wyliczanie wartosci hbc dla poszczegolnych scian
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					elem.surface[i][j][k] = detjCnt(net, i) * (pom1[i].PC1[j][k] + pom1[i].PC2[j][k] + pom1[i].PC3[j][k]);
				}
			}
		}
	}
}

void cmatrixCnt(grid& net, element4_2D elem, double ro, double cp) {
	double (*pFunc[4])(double, double);//Funkcje ksztaltu
	pFunc[0] = N1;
	pFunc[1] = N2;
	pFunc[2] = N3;
	pFunc[3] = N4;

	vec2D arr1;
	arr1.resize(elem.pointsNumber);
	for (int i = 0; i < elem.pointsNumber; ++i) {
		arr1[i].resize(4);
	}

	//wyliczanie wartosci funkcji lokalnych dla punktow calkowania
	for (int i = 0; i < elem.pointsNumber; ++i) {
		for (int j = 0; j < 4; ++j) {
			arr1[i][j] = pFunc[j](elem.wKsi[i], elem.wEta[i]);// *elem.ak[i % elem.N];
		}
	}

	//Wyliczanie macierzy C dla elementu skonczonego
	for (int i = 0; i < net.nE; ++i) {
		for (int j = 0; j < elem.pointsNumber; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < 4; ++l) {
					net.elements[i].Cmatrix[k][l] += ro * cp * elem.ak2[j] * (arr1[j][k] * arr1[j][l]) * net.elements[i].jak.det[j];
				}
			}
		}
	}
}

element4_2D derivate(int N) {
	element4_2D elem = { N };
	double (*pFunc[4])(double);//Pochodne po dEta
	pFunc[0] = dN1;
	pFunc[1] = dN2;
	pFunc[2] = dN3;
	pFunc[3] = dN4;

	double (*pFunc2[4])(double);//Pochodne po dKsi
	pFunc2[0] = dN1;
	pFunc2[1] = dN4;
	pFunc2[2] = dN3;
	pFunc2[3] = dN2;

	for (int i = 0; i < elem.pointsNumber; ++i) {
		for (int j = 0; j < 4; ++j) {
			elem.dKsi[i][j] = pFunc[j](elem.wEta[i]);
			elem.dEta[i][j] = pFunc2[j](elem.wKsi[i]);
		}
	}

	//elem.showDerivates();

	std::cout << std::endl;

	return elem;
}
#endif	
