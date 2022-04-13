function [KK_phase,F_phase]=Global_stiffness_thermal_resi_M(W,EL,phase_n_i,T_n_i,T_n,infoPHF)

%%===================================
E=infoPHF.E;
numelem=infoPHF.numelem;
element=infoPHF.element;
numnode=infoPHF.numnode;
Ck=infoPHF.Ck;
eta=infoPHF.eta;
rho_C_p=infoPHF.rho_C_p;
tau=infoPHF.tau;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
PI_vec = zeros(numelem*nnz_KE,1);
PJ_vec = PI_vec;
K_vec = zeros(numelem*nnz_KE,1);
F_phase = zeros(numnode,1);
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele,1).sctr;          % element connectivity
    nn = length(sctr) ;         % number of nodes per element
    phase_e=phase_n_i(sctr); % element phase at the respective nodes [s1 s2 s3 s4]
    T_e=T_n_i(sctr);
    % at time n
    T_e_n=T_n(sctr);
    % -------------------------------------------------------
    K_p_e=zeros(nn,nn);
    f_int_p=zeros(length(EL(ele,1).sctr'),1);
    % ---------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        B_phase=EL(ele,1).GaussValues(kk,1).B_phase;
        B_T=B_phase;
        %--------------------------------------
        T_h= N'*T_e;
        T_h_n= N'*T_e_n;
        d_h=N' * phase_e;
        Grad_T=B_T*T_e;
        g_c=(1-d_h)^(2); 
        %--------------------------------------
        % calculating Stress // Elastic Modulus:
        %--------------------------------------
        [CC,H_T]  = get_elasticity_tensor_thet(g_c,Grad_T,Ck); 
        %--------------------------------------
        % calculating residual//tangent
        %--------------------------------------
        f_micro_force1=-B_T'*(tau*H_T)+N*(rho_C_p*(T_h-T_h_n));   
        f_int_p=f_int_p+ (f_micro_force1*W(kk)*detJ);
        %--------------------------------------    
        K_phase=-B_T'*(tau*CC)*B_T+N*(rho_C_p)*N';
        K_p_e=K_p_e+(K_phase*W(kk)*detJ);
        %--------------------------------------
    end  % end of looping on GPs
    % Assemble
    PI_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = repmat(sctr',num_DOF,1);
    tmp = repmat(sctr',1,num_DOF)';
    PJ_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = tmp(:);
    K_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = K_p_e(:);
    F_phase(sctr,1)= F_phase(sctr,1) + f_int_p;
    %-----------------------------------------------------------
end % end of looping on elements
KK_phase = sparse(PI_vec,PJ_vec,K_vec,numnode,numnode);
end
%=========================================================================
function [CC,H_T]  = get_elasticity_tensor_thet(g_d,Grad_T,K_0)

K_cond=g_d*K_0;
K_tot=K_cond;
CC=-K_tot;
H_T=-K_tot*Grad_T;

end









