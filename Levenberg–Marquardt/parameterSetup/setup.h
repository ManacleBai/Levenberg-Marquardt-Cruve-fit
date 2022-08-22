#pragma once


/*parameter setup*/
#define GuessNumber 4
#define LoopMax 100
#define Lamda 0.1
#define Precision 100
#define data_length 44 


/*target function*/
#define F(x) P[0]+(P[1]/(1+exp((P[2]-OBS_X[i])/P[3])))


/*Jacobian Matrix*/
#define rep (exp((G[2] - OBS_X[i]) / G[3]))
#define Dfa 1
#define Dfb 1 / (rep + 1)
#define Dfc -(G[1] * rep) / (G[3] * (pow(rep + 1, 2)))
#define Dfd (G[1] * rep * (G[2] - OBS_X[i])) / ((G[3] * G[3]) * (pow(rep + 1, 2)))






#define nextline  cout << endl