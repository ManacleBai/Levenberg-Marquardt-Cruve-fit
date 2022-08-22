
#include "LM.h"
#include <time.h>
#include <stdlib.h>
#include "math.h"
#include <iostream>





void LM::LMFIT()
{
    /*initialization*/
    uint8_t number = 0;
    AllocateMatric(&J, data_L, parameter, number++);
    AllocateMatric(&JT, parameter, data_L, number++);
    AllocateMatric(&H, parameter, parameter, number++);
    AllocateMatric(&H_lm, parameter, parameter, number++);
    AllocateMatric(&TempArray, parameter, parameter, number++);
    /*end*/

    /*parameter setup*/
    uint16_t it = 0; // iterial max setup

    while (it < loopmax)
    {
        if (upJmatrix)
        {
            CaluateJMatric();
            fun(G);
            residual(); // calculate r(x) = y - yi;
            Hessian();  // H = J^t * J ; H = H + lamda;

            if (it == 0) VectorDot(residualArray, residualArray, &pre_e,data_L); // calculate sqrt(sum(y-yi));

        }
        else
        {
            fun(G);
            residual();
        }


        CaluateStep(); 
        /*restart calculate residual by new parameter*/
        upGuess();
        fun(temp_G);
        residual();
        VectorDot(residualArray, residualArray, &new_e, data_L);
        /*end*/

        /*calculate delta => new parameter - old parameter
        * delta : check new parameter is different from old parameter
        */
        CaluateDelta();
        /*end*/

        upJmatrix = optimizer();
        if (endLoop) break;
        it++;
    }
    cout << "********************************************************************\n";
    cout << "Iterate time : " << it << endl;
    cout << "********************************************************************\n" << endl;
    cout << "********************************************************************\n";
    cout << "Fitting parameter : ";
    MatrixPrint(&temp_G, 1, parameter);
    cout << "********************************************************************\n";
}

void LM::Init(uint16_t _loopmax, double _lamda, double _precision,uint16_t _data_L, uint8_t _parameter)
{
    loopmax = _loopmax;
    lamda = _lamda;
    data_L = _data_L;
    parameter = _parameter;
    upJmatrix = 1;
    precision = _precision;


    this->endLoop = 0;
    this->pre_e = 0;
    this->new_e = 0;
    this->delta = 0;
    this->predicetY = (double*)malloc(data_L * sizeof(double));
    this->residualArray = (double*)malloc(data_L * sizeof(double));
    this->step = (double*)malloc(_parameter * sizeof(double));
    this->temp_G = (double*)malloc(_parameter * sizeof(double));
}
void LM::DataInput(double* X, double* Y)
{
    OBS_X = (double*)malloc(data_L*sizeof(double));
    OBS_Y = (double*)malloc(data_L * sizeof(double));
    if (OBS_X != nullptr && OBS_Y != nullptr)
    {
        memcpy(OBS_X, X, data_L * sizeof(double));
        memcpy(OBS_Y, Y, data_L * sizeof(double));
    }
    else
    {
        cout << "failure";
    }
}
void LM::Guessparemter(double* _G)
{
    G = (double*)malloc(parameter * sizeof(double));
    if (G != nullptr)
    {
        memcpy(G, _G, parameter * sizeof(double));
    }

}
void LM::AllocateMatric(double*** _array, int row, int col,uint8_t number)
{
    /*
    *_array = new double*[row];
    for (int i = 0; i < row; ++i) {
        (*_array)[i] = new double[col];
    }
    */
    *_array = (double**)malloc(row * sizeof(double*));
    for (int i = 0; i < row; i++) {
        (*_array)[i] = (double*)malloc(col * sizeof(double));
    }

    switch (number)
    {
    case 0:
        this->J = J;
        for (uint16_t i = 0; i < row; i++)
        {
            for (uint16_t j = 0; j < col; j++)
                J[i][j] = 0;
        }
        break;
    case 1:
        this->JT = JT;
        for (uint16_t i = 0; i < row; i++)
        {
            for (uint16_t j = 0; j < col; j++)
                JT[i][j] = 0;
        }
        break;
    case 2:
        this->H = H;
        for (uint16_t i = 0; i < row; i++)
        {
            for (uint16_t j = 0; j < col; j++)
                H[i][j] = 0;
        }
        break;
    case 3:
        this->H_lm = H_lm;
        for (uint16_t i = 0; i < row; i++)
        {
            for (uint16_t j = 0; j < col; j++)
                H_lm[i][j] = 0;
        }
        break;
    case 4:
        this->TempArray = TempArray;
        for (uint16_t i = 0; i < row; i++)
        {
            for (uint16_t j = 0; j < col; j++)
                TempArray[i][j] = 0;
        }
        break;
    }

};
void LM::CaluateJMatric()
{
    /*calculate J maxtrix
    partial derivative for function
    exp :
          f(x) = ax^2 + bx + c
          dev(x,a);dev(x,b);dev(x,c)
    so, J = [ x^2, x, 1 ]
    */

    // calculate J maxtrix ( target function is 'a+(b/[1+exp(-(c-x)/k)])' [four parameter]

    for (uint16_t i = 0; i < data_L; i++)
    {
        switch (parameter)
        {
        case 1:
            J[i][0] = Dfa;
        case 2:
            J[i][0] = Dfa;
            J[i][1] = Dfb;
        case 3:
            J[i][0] = Dfa;
            J[i][1] = Dfb;
            J[i][2] = Dfc;
        case 4:
            J[i][0] = Dfa;
            J[i][1] = Dfb;
            J[i][2] = Dfc;
            J[i][3] = Dfd;
            break;
        default:
            cout << "failure";
            break;
        }

        //J[i][2] = -(G[1] * exp((G[2] - OBS_X[i]) / G[3]))/ (G[3] *(pow((exp((G[2] - OBS_X[i])/G[3])),2)));
        //J[i][3] = ((G[1] * exp((G[2] - OBS_X[i]) / G[3]))*(G[2] * OBS_X[i]))/((G[3]*G[3])* (pow(G[3] * (exp((G[2] - OBS_X[i]))), 2)));
    }

};
void LM::Hessian()
{
    this->JT = Transpose(J,JT, data_L, parameter); //J^T
    this->H = Multiple(JT,J,H, parameter, data_L, parameter); // H = J^t * J;
};
void LM::fun(double *P)
{
    /*target function is 'a+(b/[1+exp(-(c-x)/k)])' [four parameter]*/

    for (uint16_t i = 0; i < data_L; i++)
    {
        predicetY[i] = F(x);
    }


}
void LM::residual()
{
    for (uint16_t i = 0; i < data_L; i++)
    {
        residualArray[i] = OBS_Y[i] - predicetY[i];
    }
}
void LM::CaluateStep()
{

    double** eye = UnitMartix(parameter, parameter, lamda);

    MatrixAdd(H, H_lm,eye, parameter, parameter);
    H_lm = InverseMatrix(H_lm,parameter);
    Multiple(JT, &residualArray, TempArray, parameter, data_L, 1);
    step = *Multiple(H_lm, TempArray, &step, parameter, parameter,1);
    MatrixPrint(&step, 1, parameter);
    spaceFree(eye, parameter);
}
void LM::upGuess()
{
    memcpy(temp_G,G,parameter*sizeof(double));
    for (int i = 0; i < parameter; i++)
    {
        temp_G[i] += step[i];
    }
    this->temp_G = temp_G;
}
void LM::CaluateDelta()
{
    delta = 0;
    int value = 0;
    for (int i = 0; i < parameter; i++)
    {
        value = temp_G[i] - G[i];
        if (value < 0) value *= -1;

        delta += value;
    }

}
uint8_t LM::optimizer()
{
    if (new_e < pre_e)
    {
        if (new_e < precision || delta < 0.01) endLoop = 1;
        else
        {
            lamda /= 3;
            pre_e = new_e;
            memcpy(G, temp_G, parameter * sizeof(double));
        }
        return 1;
    }
    else
    {
        lamda *= 3;
        return 0;
    }

    return 0;
}
void LM::upThis()
{
    this->J = J;
    this->H = H;
    this->JT = JT;
}
void LM::spaceFree(double** p, int row)
{
    for (int i = 0; i < row; ++i) {
        delete p[i];
    }
    delete[] p;
};
