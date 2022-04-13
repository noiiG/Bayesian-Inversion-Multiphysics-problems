%=========================================
% @copyright(2020)
% This code is for PHF
% ** Dr. Ing Nima Noii **
% Institute of Contiunum Mechancis 
% (IKM-LUH, Hannover, DE)
%=========================================
function [Plot_U,Plot_P]=main(E,Ndof,G_f,ls,lenB,nel,gamma_c,k)
info_preprocess;
info_BC;
info_load;
initialaization_null;
info_matrirx=info_material(E,Ndof,G_f,ls,gamma_c,k);
initialaization;
%=========================================
[EL,gamma]=Global_resi_Initial(Q, W,elemType,infoPHF);
%=========================================
%tic;
for step = 1 :nstep  % new time steps
initialaization_new_step;
    for iter = 1 : 3
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
    %disp(['  S_i: ', num2str(iter)]);
    if flag_terminate==1
       break;
    end
    gamma=cal_fat_field(W,EL,u_n1_i,u_n1,phase_n1_i,phase_n1,gamma,infoPHF);
    save_results;
    %=========================================
end % end of load steps
%toc;
%=========================================
clear x
load('xxx.mat')


Plot_P=Plot_P(x);
end
