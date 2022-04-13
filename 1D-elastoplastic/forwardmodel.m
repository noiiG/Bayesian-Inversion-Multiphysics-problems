%=========================================
% Private FEM code for 2D nonlinear problems
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ H       ] hardening modulus ......................[MPa]
% [ 3] [ sigmaY  ] yield strength .........................[MPa]

%=========================================
function [Plot_P]=forwardmodel(x)
%-------------------------
E=x(1);
H  = x(2);       
sigmaY=x(3); 
%-------------------------
lenB=1;
nel=100;
h=lenB/nel;
Ndof=1;
%=========================================

info_preprocess;
info_BC;
info_load;
initialaization_null;
info_matrirx=info_material(E,H,sigmaY,Ndof);
initialaization;
flag_terminate=0;
%=========================================
[EL,KK_n1_i,Eplas_n]=Global_stiffness_Initial_cal(Q, W,elemType,I_vec, J_vec,infoPHF);
%=========================================
tic;
for step = 1 :nstep  % new time steps
initialaization_new_step;    
    %=========================================
    % 1. solving the momentum equation
    %-----------------------------------------
    elasticity_PDE;
    save_results;
    %=========================================
end % end of load steps
toc;
%=========================================
