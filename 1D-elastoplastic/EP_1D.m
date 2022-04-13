%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% Accepted for publication in Archives of Computational Methods in Engineering

% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, sand Peter Wriggers

% This script is related to 1D_elastoplasticity problem.

% @copyright(2022) LUH Hannover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ H       ] hardening modulus ......................[MPa]
% [ 3] [ sigmaY  ] yield strength .........................[MPa]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% To approximate the material properties:
% [E,H,sigmaY]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('reference.mat') 

% range of the parameters based on the prior density
pmin=[50000  200 200];
pmax=[100000 400 400];
inp.range=[pmin' pmax'];   

% number of iteration in MCMC                                
inp.nsamples=1000;          

% initial covariance                                
inp.icov=[1e7 200 200]'.*eye(size(inp.range,1));           

% initial guesss of the parameters based on the prior
inp.theta0=[60000 200 300];

% std 
inp.sigma=1;    

% measurement/ reference observation                                
inp.measurement=reference;  

% The starting point of the Kalman MCMC                                
inp.Kalmans=10;      

% assumed measurement error for Kalman MCMC                                      
inp.me=1e-1;                                               


% The DRAM algorithm 
[results] = DRAM(inp);
fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))

% The EnKF algorithm 
[results] = EnKF(inp);
fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))
