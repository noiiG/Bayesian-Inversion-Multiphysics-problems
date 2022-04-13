%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% Accepted for publication in Archives of Computational Methods in Engineering

% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, and Peter Wriggers

% This script is related to 1D_heat.

% @copyright(2022) LUH Hannover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material and model parameters

% [ 2] [ R ]       Heat capacity x density ........................[MPa/K]
% [ 3] [ lambda      ] Thermal conductivity .....................[N/K s]
% [ 4] [ v      ] Thermal velocity .....................[N/K s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% To approximate the material properties:
% [lambda,v]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
load('reference.mat')

% range of the parameters based on the prior density
pmin=[0.5 5];
pmax=[2 20];
inp.range=[pmin' pmax'];                       

% number of iteration in MCMC
inp.nsamples=1000;                             

% initial covariance 
inp.icov=[0.1 1].*eye(size(inp.range,1));  

% initial guesss of the parameters based on the prior    
inp.theta0=[1 10];                             

% std
inp.sigma=0.05;           

% measurement/ reference observation                    
inp.measurement=reference;   

% The starting point of the Kalman MCMC                  
inp.Kalmans=100;                       

% assumed measurement error for Kalman MCMC        
inp.me=1e-1;                                  

% The EnKF technique
[results] = EnKF(inp);
fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))

 
