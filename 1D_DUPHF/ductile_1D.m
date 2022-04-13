%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% Accepted for publication in Archives of Computational Methods in Engineering

% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, and Peter Wriggers

% This script is related to 1D_ductile fracture.

% @copyright(2022) LUH Hannover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ H       ] initial hardening modulus ..............[MPa]
% [ 3] [ sigmaY  ] initial yield strength .................[MPa]
% [ 4] [ l_p     ] palstic length scale ...................[mm]
% [ 5] [ psi_c   ] critical fracture energy ...............[MPa]
% [ 6] [ l_f     ] fracture length scale ..................[mm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% To approximate the material properties:
% [E,sigmaY,H,psi_c]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load('reference.mat')

% range of the parameters based on the prior density
pmin=[30000 300 200 10];
pmax=[100000 400 400 50];
inp.range=[pmin' pmax'];

% number of iteration in MCMC                                
inp.nsamples=1000;           

% initial covariance                            
inp.icov=[1e7 100 100 10]'.*eye(size(inp.range,1));  

% initial guesss of the parameters based on the prior   
inp.theta0=[50000 300 300 20];                         

% std
inp.sigma=0.1;

% measurement/ reference observation                                
inp.measurement=reference;          

% The starting point of the Kalman MCMC                    
inp.Kalmans=100;       

% assumed measurement error for Kalman MCMC                                 
inp.me=1e-1;                                            


% The EnKF technique
[results] = EnKF(inp);
fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))
