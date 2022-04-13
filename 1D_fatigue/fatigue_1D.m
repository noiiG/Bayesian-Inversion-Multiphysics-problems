%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% Accepted for publication in Archives of Computational Methods in Engineering

% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, and Peter Wriggers

% This script is related to 1D_fatigue problem.

% @copyright(2022) LUH Hannover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ G_f     ] fracture toughness .....................[MPa mm]
% [ 3] [ l_f     ] fracture length scale ..................[mm]
% [ 4] [ gamma_c ] fatigue threshold ......................[MPa]
% [ 5] [ k       ] fatigue degradation parameter ..........[-]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% To approximate the material properties:
% [E,k,G_f,gamma_c]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('reference.mat')
load('xxx.mat')

% range of the parameters based on the prior density
pmin=[50000  0.2 0.02 1e-3];
pmax=[100000 0.6 0.04 1e-2];
inp.range=[pmin' pmax'];   

% number of iteration in MCMC                                       
inp.nsamples=1000;       
                                         
% initial covariance 
inp.icov=[1e7 1e-3 1e-4 1e-6]'.*eye(size(inp.range,1));          

% initial guesss of the parameters based on the prior
inp.theta0=[60000 0.3 0.03 5e-3];                    

% measurement/ reference observation             
inp.sigma=0.1;                      

% The starting point of the Kalman MCMC           
inp.measurement=reference(x);  

inp.Kalmans=10;    

% assumed measurement error for Kalman MCMC                                             
inp.me=1e-1;                                                    

% The inferred parameters using DRAM algorithm
[results] = DRAM(inp);
fprintf('The median of the posterior is:%d\n',median((results.MCMC)'))
  

