#pragma once

//void MatrixPrint(double** matrix, int row, int col);
class MatrixOperation
{
public:
	void MatrixPrint(double** matrix, int i_size, int j_size);
	double** Transpose(double** matrix,double **buffer, int row, int col);
	double** Multiple(double** A, double **B,double **buffer,int A_row,int B_row,int B_col);
	void VectorDot(double* A, double* B, double *buffer,int length);
	double** UnitMartix(int row,int col,double value);
	void MatrixAdd(double** augend, double** A, double** B, int row, int col);
	double** InverseMatrix(double** A, int element);
};