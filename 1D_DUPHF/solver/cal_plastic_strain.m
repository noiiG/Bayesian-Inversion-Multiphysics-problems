function [Eplas_new]=cal_plastic_strain(W,EL,d_n,phase_n1,...
    alpha_n_i,alpha_n,Eplas,infoPHF)

%===================================
Eplas_new=Eplas;
%===================================
element=infoPHF.element;
E=infoPHF.E;
numelem=infoPHF.numelem;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele,1).sctr;          % element connectivity
    nn = length(sctr) ;         % number of nodes per element
    disp_u=d_n(sctr,:);    % element disp at the respective nodes[u1 u2 u3 u4 u5 u6 u7 u8]
    phase=phase_n1(sctr); % element phase at the respective nodes [s1 s2 s3 s4]
    alphapn=alpha_n(sctr);
    alphap=alpha_n_i(sctr);
    % -------------------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        B_phase=EL(ele,1).GaussValues(kk,1).B_phase;
        %--------------------------------------
        strain_Pn=Eplas(ele,1).data(kk,1).s;
        % calculating field values at GP at current iteration
        strain=Bfem*disp_u;
        %--------------------------------------
        [normS]= cal_norm_J2S(strain,strain_Pn);
        del_alpha=N'*(alphap-alphapn);
        strain_p=strain_Pn+del_alpha*normS;
        Eplas_new(ele,1).data(kk,1).s=strain_p;
        %--------------------------------------
    end  % end of looping on GPs
end % end of looping on elements
end
%=========================================================================
function[snorm]= cal_norm_J2S(strain,strain_Pn)
%------------------------
snorm=sign(strain-strain_Pn);
end
%=========================================================================

