function [Ce,C_vol,C_dev]=linear_elasticity(lamda,mu,K)

Ce = [lamda+(2*mu)   lamda           0;
    lamda          lamda+(2*mu)    0;
    0              0            mu];

C_vol =K*[1   1    0;
    1   1    0;
    0   0    0];

C_dev =[2/3  -1/3  0;
    -1/3   2/3  0;
    0     0   1/2];
