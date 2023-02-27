#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <cmath>

#include "headers\matrix.hpp"

// Количество узлов
int countPoints;
// Количество треугольников 
int countTriangles;
// Количество граничных линий 
int countLines; 

// Коэффициенты фильтрации
double Kxx = 40.0; 
double Kyy = 20.0;

// Просачивание воды (вдоль реки)
double q = -108.0; 
// Мощности насоса
double P1 = -1200, P2 = -2400; 
// Расположение координат насосов
double coord_x1 = 25.0, coord_y1 = 10.0, coord_x2 = 8.0, coord_y2 = 5.0;
// Напор
const double fi = 200.0;

// Координаты узлов
std::vector<std::vector<double>> points;
// Номера узлов треугольника
std::vector<std::vector<int>> triangles;
// Номер границы и узлы границы
std::vector<std::vector <int>> lines;
// Площадь треугольника
std::vector<double> squareTriangles;
	
// Длина линии 
double lineLength; // Длина линии

void setPoints(std::ifstream &file) {    
    points.resize(countPoints);
    for (auto i = 0; i < countPoints; i++) {
        points[i].resize(3);
		file >> points[i][0];
		file >> points[i][1];
		file >> points[i][2];
    }
}

void setTriangles(std::ifstream &file) {
    triangles.resize(countTriangles);
    for (auto i = 0; i < countTriangles; i++) {
		triangles[i].resize(3);
		// Записываем ID материала
		file >> triangles[i][0];
		// Перезаписываем ID материала на ID вершины треугольника
		file >> triangles[i][0]; 
		file >> triangles[i][1];
		file >> triangles[i][2];
	}
}

void setLines(std::ifstream &file) {
    lines.resize(countLines);
	for (int i = 0; i < countLines; i++) {
		lines[i].resize(3);
		file >> lines[i][0]; // Номер границы
		file >> lines[i][1];
		file >> lines[i][2];
	}
}

Matrix setCoordinatesForGlobalMatrix(Matrix C, Matrix D, Matrix K) {
    for (int k = 0; k < countTriangles; k++) {
		// Координаты КЭ (первый столбец из 1)
        Matrix C(3,3);
        C.SetCoordinatesFE(k, triangles, points);
        // C.Show("Матрица C");

		squareTriangles[k] = 1 / 2.0 * fabs(C.Det(3));

		Matrix inv_C = C.Inversion();
        // inv_C.Show("Обратная матрица С");

        Matrix B(2,3);
		// Хранит коэффициенты функции формы
        B.SetBMatrix(inv_C);
        // B.Show("Матрица B");
  

        Matrix B_tr = B.Transpose();
        // B_tr.Show("Транспонированная матрица B");
        
		// Умножение матриц B^T*D*B 
        Matrix BDB = B_tr.MultiplicationByMatrix(D.MultiplicationByMatrix(B));
		BDB.MultiplicationByNumber(squareTriangles[k]);
        // BDB.Show("Матрица B^T*D*B");

		// Переносим в глобальную матрицу K
		for (int t = 0; t < 3; t++) {
			for (int g = 0; g < 3; g++) {
				int ind1, ind2;
				ind1 = triangles[k][t];
				ind2 = triangles[k][g];
                K.SetElem(ind1 - 1, ind2 - 1, K.GetElem(ind1 - 1, ind2 - 1) + BDB.GetElem(t, g));
			}
		}
	}
    return K;
}

double getLength(std::vector<std::vector<int>> lines, std::vector<std::vector<double>> points, int k) {
	int id1 = lines[k][1];
	int	id2 = lines[k][2];
	double length = 0.0;
	for (uint8_t i = 0; i < 2; ++i) {
		length += std::pow((points[id1 - 1][i] - points[id2 - 1][i]), 2);
	}
	return sqrt(length);
}
std::vector<double> getFuncionOfForm(double x, double y, Matrix inv_C) {
	return {inv_C.GetElem(0,0) + inv_C.GetElem(1,0) * x + inv_C.GetElem(2,0) * y, 
            inv_C.GetElem(0,1) + inv_C.GetElem(1,1) * x + inv_C.GetElem(2,1) * y, 
            inv_C.GetElem(0,2) + inv_C.GetElem(1,2) * x + inv_C.GetElem(2,2) * y};
}

int main() {
    std::string meshFile = "mesh.neu";

    std::ifstream file;
    file.open(meshFile);
	if(!file.is_open()) {
		std::cout << "Ошибка открытия файла" << std::endl;
		return 0;
	}

    file >> countPoints;
    std::cout << "Число точек = " << countPoints << std::endl;
    setPoints(file);

    // Глобальная матрица K
    Matrix K(countPoints, countPoints);

    // Глобальный вектор (правая часть)
    std::vector<double> F;
    F.resize(countPoints);

    // Глобальный вектор значений (решение)
    std::vector<double> Fi;
    Fi.resize(countPoints);

    file >> countTriangles;
    std::cout << "Число треугольников = " << countTriangles << std::endl;
    setTriangles(file);


    file >> countLines;
	std::cout << "Число граничных линий: " << countLines << std::endl;
	lines.resize(countLines);
    setLines(file);

    squareTriangles.resize(countTriangles);

    file.close();

    // Матрица свойств
    Matrix D(2,2);
    D.SetProperties(Kxx, Kyy);
    // D.Show("Матрица свойств D");

    // Матрица координат КЭ
    Matrix C(3, 3);
    
    K = setCoordinatesForGlobalMatrix(C, D, K);
    K.Show("Глобальная матрица K");

    // Граничные условия
    std::vector<int> id_bc;
	for (int k = 0; k < countLines; k++) {
        // ГУ 1 рода (задано значение fi)
		if (lines[k][0] == 1 || lines[k][0] == 3) {
			int ind1 = lines[k][1];
			int ind2 = lines[k][2];
			Fi[ind1-1] = fi;
			Fi[ind2-1] = fi;
			id_bc.push_back(ind1 - 1);
			id_bc.push_back(ind2 - 1);
		}

        // ГУ 2 рода (задано q)
		if (lines[k][0] == 5) {
			lineLength = getLength(lines, points, k);
			// Перенос в глобальный вектор нагрузки F
			for (int i = 1; i < 3; i++) {
				int ind = lines[k][i]; 
				F[ind - 1] += -(q * lineLength)/2.0;
			}
		}

	}
    auto id_boundary_condition = std::unique(id_bc.begin(), id_bc.end());
	id_bc.erase(id_boundary_condition, id_bc.end());
	Matrix C1(3,3);
	for (int k = 0; k < countTriangles; k++) {
		//Проверяем принадлежит ли насос данному КЭ
		//Получаем ID узлов для данного КЭ
		int id1 = triangles[k][0], id2 = triangles[k][1], id3 = triangles[k][2];
		// Координаты полученных узлов
		double x1 = points[id1 - 1][0], x2 = points[id2 - 1][0], x3 = points[id3 - 1][0];
		double y1 = points[id1 - 1][1], y2 = points[id2 - 1][1], y3 = points[id3 - 1][1];

		//для 1 насоса
		if (!((x1 - coord_x1) * (y2 - y1) - (x2 - x1) * (y1 - coord_y1) < 0 ||
				(x2 - coord_x1) * (y3 - y2) - (x3 - x2) * (y2 - coord_y1) < 0 ||
				(x3 - coord_x1) * (y1 - y3) - (x1 - x3) * (y3 - coord_y1) < 0)) {
			
            double square = squareTriangles[k];

			int s = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (j == 0) {
                        C1.SetElem(i,j,1);
						s++;
					}
					else {
						int id = triangles[k][s - 1] - 1;
						C1.SetElem(i,j,points[id][j - 1]);
					}
				}
			}
			std::vector<double> N = getFuncionOfForm(coord_x1, coord_y1,C1.Inversion());
			F[id1 - 1] += N[0] * P1 * square / 3.0;
			F[id2 - 1] += N[1] * P1 * square / 3.0;
			F[id3 - 1] += N[2] * P1 * square / 3.0;
		}

		//для 2 насоса
		if (!((x1 - coord_x2) * (y2 - y1) - (x2 - x1) * (y1 - coord_y2) < 0 ||
			(x2 - coord_x2) * (y3 - y2) - (x3 - x2) * (y2 - coord_y2) < 0 ||
			(x3 - coord_x2) * (y1 - y3) - (x1 - x3) * (y3 - coord_y2) < 0)) {
			double square = squareTriangles[k];

			int s = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (j == 0) {
						C1.SetElem(i, j , 1);
						s++;
					}
					else {
						int id = triangles[k][s - 1] - 1;
						C1.SetElem(i, j, points[id][j - 1]);
					}
				}
			}
			std::vector<double> N = getFuncionOfForm(coord_x2, coord_y2, C1.Inversion());
			F[id1 - 1] += N[0] * P2 * square / 3.0;
			F[id2 - 1] += N[1] * P2 * square / 3.0;
			F[id3 - 1] += N[2] * P2 * square / 3.0;
		}
	}

	// Преобразование системы с учетом ГУ 1 рода
	for (auto &i: id_bc) {
		F[i] = K.GetElem(i,i) * Fi[i];

		for (int j = 0; j < countPoints; j++) {
			if (j == i)
				continue;

			F[j] -= K.GetElem(i,j) * Fi[i];
			K.SetElem(i, j, 0.0);
			K.SetElem(j, i, 0.0);
		}
	}

	//Решение СЛАУ K*Fi=F
	std::cout << "Решение:" << std::endl;
	Matrix K_inv = K.Inversion();
	for (int i = 0; i < countPoints; i++) {
		for (int j = 0; j < countPoints; j++)
			if(std::find(id_bc.begin(),id_bc.end(),i)==id_bc.end())
				Fi[i] += K_inv.GetElem(i,j) * F[j];
		std::cout << Fi[i] << std::endl;
	}
	return 0;
}
