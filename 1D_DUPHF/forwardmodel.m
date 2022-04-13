%========================================= 
function [Plot_P]=forwardmodel(x)
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ H       ] initial hardening modulus ..............[MPa]
% [ 3] [ sigmaY  ] initial yield strength .................[MPa]
% [ 4] [ l_p     ] palstic length scale ...................[mm]
% [ 5] [ psi_c   ] critical fracture energy ...............[MPa]
% [ 6] [ l_f     ] fracture length scale ..................[mm]

%=========================================
E=x(1);
sigmaY=x(2);
H=x(3);
psi_c=x(4);
%--------------------------
l_p=10/sqrt(sigmaY); %2*ls
ls=0.55/(sqrt(2)*sqrt(psi_c));%2*h;         
%--------------------------
EP_bar_crit=0.1;
Ndof=1;
p_stable=0;
chi=55;
G_f=0;
%--------------------------
theta=45;
lenB=1;
nel=100;
h=lenB/nel;
at=1;
%--------------------------
irr_phase=1;      % 1 for ad hoc, 2 for penalty, 3 for H.F., 4 for D 
irr_alpha=1; 
%=========================================

info_preprocess;
info_BC;
info_load;
initialaization_null;
info_matrirx=info_material(E,psi_c,H,sigmaY,EP_bar_crit,p_stable,Ndof,chi,l_p,G_f,at,irr_phase,irr_alpha,ls,theta);
initialaization;
%=========================================
[EL,Eplas_n]=Global_resi_Initial(Q, W,elemType,infoPHF);
%=========================================
for step = 1 :nstep  % new time steps
initialaization_new_step;
flag_comp=0;
    for iter = 1 : 300
%         disp(['  S_i: ', num2str(iter)]);
        %=========================================
        % 1. solving the momentum equation
        %--------------------------------------------------
        elasticity_PDE;
        %--------------------------------------------------        
        % 2. solving the plasticity equation
        %--------------------------------------------------
        plasticity_PDE;
        %--------------------------------------------------
        % 3. solving the phase field equation
        %--------------------------------------------------  
        phase_field_PDE; 
        %--------------------------------------------------
        % Convergence
        %--------------------------------------------------
         if max([norm(u_n1_i-u_n1_0,inf),norm(alpha_n_i-alpha_n_0,inf),norm(phase_n1_i-phase_n1_0,inf)])<1e-5
             break
         else
             u_n1_0=u_n1_i;    
             alpha_n_0=alpha_n_i;
             phase_n1_0=phase_n1_i;
        end  
         %--------------------------------------------------
    end     
    %=========================================
    save_results;
    %=========================================
end % end of load steps
%=========================================

