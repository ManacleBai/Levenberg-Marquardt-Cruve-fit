// LMFIT.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//
#include <iostream>
#include <string>
#include "parameterSetup/setup.h"
#include "src/LM.h"

using namespace std;



int main()
{

    LM a;
    double testdataX[] =
    {
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, \
        21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38, \
        39,40,41,42,43
    };
    double testdataY[] =
    {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,33,62,108,173,267,395,552,731,925,\
        1125,1322,1514,1696,1862,2009,2130,2224,2298,2338,2356,2356,2356,\
        2356,2356,2370,2387,2408,2432,2457
    };

    double Guessparameter[] = { 0,2457,0.1,0.1 }; // a,b,c,d
    double START, END; START = clock();
    a.Init(LoopMax, Lamda, Precision, data_length, GuessNumber);
    a.DataInput(testdataX, testdataY);
    a.Guessparemter(Guessparameter);
    a.LMFIT();

    cout << "********************************************************************\n";
    cout << "process Time : ";
    END = clock();
    cout << (END - START) / CLOCKS_PER_SEC << "ms" << endl;
    cout << "********************************************************************\n";
    return 0;
}


