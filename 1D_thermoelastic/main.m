%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is related to the following paper:
% Bayesian inversion with open-source codes for various one-dimensional model problems in computational mechanics
% submitted to Archives of Computational Methods in Engineering

% The authors are:
% Nima Noii, Amirreza Khodadadian, Jacinto Ulloa, Fadi Aldakheel, Thomas Wick, and Peter Wriggers

% This script is related to 1D_fatigue problem.

% @copyright(2022) LUH Hannover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Plot_U,Plot_P]=main(E,Ndof,G_f,irr_phase,ls,Calpha,rho_C_p,Ck,TimeF,lenB,nel)
info_preprocess;
info_BC;
info_load;
initialaization_null;
info_matrirx=info_material(E,Ndof,G_f,irr_phase,ls,Calpha,rho_C_p,Ck);
initialaization;
%=========================================
[EL,ET_n]=Global_resi_Initial(Q, W,elemType,infoPHF);
%=========================================
%tic;
for step = 1 :nstep  % new time steps
initialaization_new_step;
    for iter = 1 : 2
%         disp(['  S_i: ', num2str(iter)]);
        %=========================================
        % 1. solving the momentum equation
        %--------------------------------------------------
        elasticity_PDE;
        %--------------------------------------------------        
        % 2. solving the plasticity equation
        %--------------------------------------------------
        thermal_PDE;
        %--------------------------------------------------
%         % Convergence
%         %--------------------------------------------------
%         if max([norm(u_n1_i-u_n1_0,inf),norm(T_n_i-T_n_0,inf),norm(phase_n1_i-phase_n1_0,inf)])<0.0001
%             break
%         else
%             u_n1_0=u_n1_i;    
%             T_n_0=T_n_i;
%             phase_n1_0=phase_n1_i;
%         end  
        %--------------------------------------------------
    end  
    %--------------------------------------------------
    % 3. solving the phase field equation
    %--------------------------------------------------
    phase_field_PDE; 
    %--------------------------------------------------
    %=========================================
    %disp(['  S_i: ', num2str(iter)]);
    if flag_terminate==1
       break;
    end
    save_results;
    %=========================================
end % end of load steps
%toc;
%=========================================
%save('Result.mat');
%PLOT_PHF;
%=========================================
