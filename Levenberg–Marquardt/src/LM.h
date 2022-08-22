#pragma once

#include "MatrixOperation.h"
#include "/Users/88698/Desktop/Code/Levenberg¡VMarquardt/parameterSetup/setup.h"

using namespace std;



typedef  unsigned short int uint16_t;
typedef unsigned char uint8_t;


class LM : public MatrixOperation
{
private:
	uint16_t data_L;
	uint8_t parameter;
	uint16_t loopmax;
	double lamda;
	uint8_t upJmatrix;
	uint8_t endLoop;
	double* G; //Guessparameter
	double* temp_G; //new_guessparaemter
	double* OBS_X;
	double* OBS_Y;
	double* predicetY;
	double* step;
	double new_e; //this time calculate(simulation) sum(residual)
	double precision;
	double delta;
protected :
	double pre_e; //last time calculate sum(residual)
public :
	double** J;
	double** JT;
	double** H;
	double** H_lm;
	double** TempArray;
	double* residualArray;
	void LMFIT();
	void Init(uint16_t _loopmax,double _lamda,double _precision,uint16_t _data_L, uint8_t _parameter);
	void DataInput(double *X, double *Y);
	void Guessparemter(double* _G);
	void AllocateMatric(double*** _array,int row,int col, uint8_t number);
	void CaluateJMatric();
	void Hessian();
	void fun(double *P);
	void residual();
	void CaluateStep();
	void upGuess();
	void CaluateDelta();
	uint8_t optimizer();
	void upThis();
	void spaceFree(double **p,int row);
};