function [model] = forwardmodel(x)
%=========================================
% 1D thermal fracture
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ G_c     ] fracture toughness .....................[MPa mm]
% [ 3] [ l_f     ] fracture length scale ..................[mm]
% [ 4] [ rho_C_p ] Heat capacity x density ........................[MPa/K]
% [ 5] [ Ck      ] Thermal conductivity .....................[N/K s]
% [ 6] [ Calpha  ] Linear thermal expansion ..................[1/K]

%=========================================
psi_c=30;
Ndof=1;
G_f=x(1);
%-----------------
E = x(2);               
ls=0.55/(sqrt(2)*sqrt(psi_c));
%-----------------
rho_C_p=0.775*2450;
%Ck=2e-2;
%Calpha=8e-4;
Ck=x(3);             
Calpha=x(4);
%-----------------
irr_phase=1;       
%-----------------
lenB=pi;
nel=99;
h=lenB/nel;
%-----------------
TimeF=1e5;
% To check: Tau>rho_C_p*h^2/(6*Ck)
%-----------------
%=========================================
[~,model]=main(E,Ndof,G_f,irr_phase,ls,Calpha,rho_C_p,Ck,TimeF,lenB,nel);
%=========================================
end

