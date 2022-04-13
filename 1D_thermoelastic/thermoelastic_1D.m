 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% Accepted for publication in Archives of Computational Methods in Engineering


% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, and Peter Wriggers

% This script is related to thermoelastic_1D fracture.

% @copyright(2022) LUH Hannover


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ G_f     ] fracture toughness .....................[MPa mm]
% [ 3] [ l_f     ] fracture length scale ..................[mm]
% [ 4] [ rho_C_p ] Heat capacity x density ........................[MPa/K]
% [ 5] [ Ck      ] Thermal conductivity .....................[N/K s]
% [ 6] [ Calpha  ] Linear thermal expansion ..................[1/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% To approximate the material properties:
% [G_f,E,Ck,Calpha]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 clc
 clear
 
load('reference.mat')
pmin=[0.01 200 1e-2 7e-4];
pmax=[0.1  400 1e-1 5e-3];
 
inp.range=[pmin' pmax'];                                       % range of the parameters based on the prior density
inp.nsamples=1000;                                             % number of iteration in MCMC
inp.icov=[0.0001 100 1e-1 1e-7].*eye(size(inp.range,1));       % initial covariance 
inp.theta0=[0.05 300 5e-2 1e-3];                               % initial guesss of the parameters based on the prior
inp.sigma=0.1;                                
inp.measurement=reference;                                     % measurement/ reference observation
inp.Kalmans=20;                                                % The starting point of the Kalman MCMC
inp.me=1e-1;                                                   % assumed measurement error for Kalman MCMC

tic
[results] = EnKF(inp);
toc

fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))
fprintf('The acceptance rate is:%d\n',mean(results.accepted))
