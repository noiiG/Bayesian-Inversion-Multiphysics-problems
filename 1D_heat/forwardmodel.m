function [model] = forwardmodel(x)

%=========================================
% The heat evolution code is taken from 
% Schirén, W. (2018). Finite element method for 1D transient convective heat transfer problems. 

%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ R ]       Heat capacity x density ........................[MPa/K]
% [ 3] [ lambda      ] Thermal conductivity .....................[N/K s]
% [ 4] [ v      ] Thermal velocity .....................[N/K s]

%=========================================
 c =1.2;  
 rho=1.6;
 lambda=x(1);
 v=x(2);
R = rho*c;
%-------------------------
noe = 100; % Nr of elements
 non = noe+1; % Nr of nodes
 not = 100; % Nr of timesteps
 T = 0.2; % Total time
 dT=T/not;
 L = 1; % Total length
 elLength = L/noe; % Element length
 Theta = 2/3; % Capital theta, 0  Theta  1;
 Pe = abs(v*elLength/2/lambda); % P´eclet number
 alphaOpt = coth(Pe) - 1/Pe; %
 gamma = v/abs(v); %

 %------------FEM----------------------------------------------------

 % -----------Build Edof-------------------------
 % Element number i goes from node i to node i+1.
 Edof = zeros(noe,3);
 for i=1:noe
 Edof(i,1)=i;
 Edof(i,2)=i;
 Edof(i,3)=i+1;
 end

 % -----------Boundary conditions and empty matrices---------------
 bc = [1 1;non 0]; % Boundary conditions, at first node temp = 1C,at last node temp = 0C
 initialTemp = 0;
 an = zeros(non,1)+initialTemp; % Initial temperatures at each node
 K = zeros(non, non); % Empty element stiffness matrix
 C = zeros(non, non); % Empty element damping matrix
 A_matrix = zeros(non, non); % Matrix to store results in

 % -----------Spatial domain---------------
 Ke = lambda/elLength*[1 -1; -1 1]; % Element stiffness
 Ktilde = v/2*[-1 1;-1 1] + alphaOpt*gamma*v/2*[1 -1;-1 1]; % Element stiffness

 Ce = rho*c*elLength/6*[2 1; 1 2] + ...
alphaOpt*gamma*rho*c*elLength/4*[-1 -1;1 1]; % Element damping


 % -----------Global matrices--------------
 K_global = assem(Edof, K, Ke);
 Ktilde_global = assem(Edof, K, Ktilde);
 C_global = assem(Edof, C, Ce);

 % -----------Time domain--------------------------------------------
 % -----------Time step zero---------------
 A_matrix(:,1) = an; % Add the initial temperatures to the ...first column of A matrix
 A_matrix(1,1) = 1; % Add the boundary condition

 A_matrix_tilde(:,1) = an;
 A_matrix_tilde(1,1) = 1;
 % -----------The development over time--------------
 
 
 
 for i=1:not
 K_star=C_global/dT+Theta*(K_global+Ktilde_global);
 f_star = (C_global/dT-(K_global + Ktilde_global) *(1-Theta))*an; % Global damping with respect to time
an1 = solveq(K_star, f_star, bc); % Temperature in each node at time n+1
A_matrix(:,i+1) = an1; % Store the temperatures in column i+1
an = an1; % Update temperatures for ...time n
end
  
 model=(A_matrix(:,100));
end

