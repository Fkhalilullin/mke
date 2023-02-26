#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

class Matrix {
    private:
        std::vector<std::vector<double>> m_matrix;

    public:
        Matrix(int N, int M);

        void MultiplicationByNumber(double number);
        Matrix MultiplicationByMatrix(Matrix m);
        void AdditionByMatrix(Matrix matrix);
        double Det(int n);
        Matrix Inversion();
        Matrix Transpose();

        void SetProperties(double Kxx, double Kyy);
        void SetCoordinatesFE(int triangleID, std::vector<std::vector<int>> triangles, std::vector<std::vector<double>> points);

        double GetElem(int i, int j);
        void SetElem(int i, int j, double elem);

        void SetBMatrix(Matrix inv_C);

        void Show();
        void Show(std::string title);

    private:
        std::vector<std::vector<double>> getMatrix(){
            return this->m_matrix;
        }
};

Matrix::Matrix(int N, int M) {
    std::vector<std::vector<double>> result;

    for(auto row = 0; row < N; row++) {
        result.push_back(std::vector<double>());
        for(auto col = 0; col < M; col++) {
            result[row].push_back(0.0);
        }
    }

    this->m_matrix = result;
}

// Умножает матрицу на число
void Matrix::MultiplicationByNumber(double number) {
    size_t rowSize = this->m_matrix.size();
    size_t colSize = this->m_matrix[0].size();

    for (auto row = 0; row < rowSize; row++)
        for (auto col = 0; col < colSize; col++)
            this->m_matrix[row][col] *= number;
}

// Выводит матрицу на экран
void Matrix::Show() {
    size_t rowSize = this->m_matrix.size();
    size_t colSize = this->m_matrix[0].size();

    for (auto row = 0; row < rowSize; row++) {
        for (auto col = 0; col < colSize; col++){
            std::cout << std::setw(5) << m_matrix[row][col] << " ";
        }
        std::cout << std::endl;
    }
}

// Вычисляет сумму двух матриц
void Matrix::AdditionByMatrix(Matrix inputMatrix) {
    size_t rowSize = this->m_matrix.size();
    size_t colSize = this->m_matrix[0].size();

     std::vector<std::vector<double>> matrix = inputMatrix.getMatrix();

    for (auto row = 0; row < rowSize; row++)
        for (auto col = 0; col < colSize; col++)
            this->m_matrix[row][col] += matrix[row][col];
}

// Заполнение матрицы свойств
void Matrix::SetProperties(double Kxx, double Kyy) {
 	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (i==0 && j==0)
				this->m_matrix[i][j] = Kxx;
			if (i == 1 && j == 1)
				this->m_matrix[i][j] = Kyy;
			if (i!=j)
				this->m_matrix[i][j] = 0;
		}
	}
}

// Выводит матрицу на экран
void Matrix::Show(std::string title) {
    size_t rowSize = this->m_matrix.size();
    size_t colSize = this->m_matrix[0].size();

    std::cout << std::endl << title << ": " << std::endl;
    for (auto row = 0; row < rowSize; row++) {
        for (auto col = 0; col < colSize; col++){
            std::cout << m_matrix[row][col] << " ";
        }
        std::cout << std::endl;
    }
}

void Matrix::SetCoordinatesFE(
    int triangleID,
    std::vector<std::vector<int>> triangles, 
    std::vector<std::vector<double>> points) {
    int s = 0;
    for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
				if (j == 0) {
					this->m_matrix[i][j] = 1;
					s++;
				}
				else {
					int id = triangles[triangleID][s - 1] - 1;
					this->m_matrix[i][j] = points[id][j - 1];
				}
		}
	}
}

// Считает определитель матрицы размером n x n
double Matrix::Det(int n) {
    if (n == 1)
		return this->m_matrix[0][0];
	else if (n == 2)
		return this->m_matrix[0][0] * this->m_matrix[1][1] - 
        this->m_matrix[0][1] * this->m_matrix[1][0];
	else {
		double d = 0;
		for (int k = 0; k < n; k++) {
            Matrix m(n-1, n-1);
			for (int i = 1; i < n; i++) {
				int t = 0;
				for (int j = 0; j < n; j++) {
					if (j == k)
						continue;
					m.m_matrix[i - 1][t] = this->m_matrix[i][j];
					t++;
				}
			}

			d += pow(-1, k + 2) * this->m_matrix[0][k] * m.Det(n - 1);
            
		}
		return d;
	}
}

Matrix Matrix::Inversion() {
    double temp;
    int N = this->m_matrix.size();

    // Единичная матрица
    Matrix E(N, N);
    Matrix A = *this;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			E.m_matrix[i][j] = 0.0;

			if (i == j)
				E.m_matrix[i][j] = 1.0;
		}
    }

	for (int k = 0; k < N; k++) {
		temp = A.m_matrix[k][k];

		for (int j = 0; j < N; j++) {
			A.m_matrix[k][j] /= temp;
			E.m_matrix[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++) {
			temp = A.m_matrix[i][k];

			for (int j = 0; j < N; j++) {
				A.m_matrix[i][j] -= A.m_matrix[k][j] * temp;
				E.m_matrix[i][j] -= E.m_matrix[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--) {
		for (int i = k - 1; i >= 0; i--) {
			temp = A.m_matrix[i][k];

			for (int j = 0; j < N; j++) {
				A.m_matrix[i][j] -= A.m_matrix[k][j] * temp;
				E.m_matrix[i][j] -= E.m_matrix[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A.m_matrix[i][j] = E.m_matrix[i][j];

    return A;
}

// Берем вторую и третью строчку обратной матрицы C
void Matrix::SetBMatrix(Matrix inv_C) {
    for (int i = 0; i < 2; i++) {
	    for (int j = 0; j < 3; j++) {
				this->m_matrix[i][j] = inv_C.m_matrix[i + 1][j];
		}
	}
}

Matrix Matrix::Transpose() {
    Matrix mat_tr(this->m_matrix[0].size(), this->m_matrix.size());

    for (int i = 0; i < this->m_matrix[0].size(); i++) {
		for (int j = 0; j < this->m_matrix.size(); j++) {
				mat_tr.m_matrix[i][j] = this->m_matrix[j][i];
		}
	}
    
    return mat_tr;
}

// Умножение матрицы на матрицу
Matrix Matrix::MultiplicationByMatrix(Matrix mat) {
    int n = int(this->m_matrix.size());
    int m = int(this->m_matrix[0].size());
    int q = int(mat.m_matrix.size());
    int p = int(mat.m_matrix[0].size());

	if (m == q) {
		Matrix res(n, p);

		for (int row = 0; row < n; row++) {
			for (int col = 0; col < p; col++) {
				for (int j = 0; j < m; j++) {
					res.m_matrix[row][col] += this->m_matrix[row][j] * mat.m_matrix[j][col];
				}
			}
		}
		return res;
	}
}

double Matrix::GetElem(int i, int j) {
    return this->m_matrix[i][j];
}

void Matrix::SetElem(int i, int j, double elem) {
    this->m_matrix[i][j] = elem;
}