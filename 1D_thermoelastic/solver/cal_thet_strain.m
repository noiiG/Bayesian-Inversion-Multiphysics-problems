function [ET_new]=cal_thet_strain(W,EL,T_n_i,infoPHF)

%===================================
element=infoPHF.element;
E=infoPHF.E;
numelem=infoPHF.numelem;
Calpha=infoPHF.Calpha;
T0=infoPHF.theta0;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele,1).sctr;          % element connectivity
    T_e=T_n_i(sctr);
    T0_e=T0(sctr);
    % -------------------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        T0_gp=N'*T0_e;
        T_gp=N'*T_e;
        strain_T=Calpha*(T_gp-T0_gp);
        ET_new(ele,1).data(kk,1).s=strain_T;
        %--------------------------------------
    end  % end of looping on GPs
end % end of looping on elements
end
%=========================================================================

