%=========================================
function [Plot_P]=forwardmodel(x)
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ G_c     ] fracture toughness .....................[MPa mm]
% [ 3] [ l_f     ] fracture length scale ..................[mm]

%=========================================

psi_c=30;

E=x(1);
G_f=x(2);

ls=0.55/(sqrt(2)*sqrt(psi_c));%2*  

Ndof=1;

lenB=1;
nel=300;
h=lenB/nel;
%=========================================
info_preprocess;
info_BC;
info_load;
initialaization_null;
info_matrirx=info_material(E,Ndof,G_f,ls);
initialaization;
%=========================================
[EL]=Global_resi_Initial(Q, W,elemType,infoPHF);
%=========================================
for step = 1 :nstep  % new time steps
initialaization_new_step;
    for iter = 1 : 20
%         disp(['  S_i: ', num2str(iter)]);
        %=========================================
        % 1. solving the momentum equation
        %--------------------------------------------------
        elasticity_PDE;
        %--------------------------------------------------        
        % 2. solving the phase field equation
        %--------------------------------------------------
        phase_field_PDE; 
        %--------------------------------------------------
        % Convergence
        %--------------------------------------------------
        if max([norm(u_n1_i-u_n1_0,inf),norm(phase_n1_i-phase_n1_0,inf)])<0.0001
            break
       else
            u_n1_0=u_n1_i;    
            phase_n1_0=phase_n1_i;
        end  
        %--------------------------------------------------
    end       
    %=========================================
    save_results;
    %=========================================
end % end of load steps
%=========================================


