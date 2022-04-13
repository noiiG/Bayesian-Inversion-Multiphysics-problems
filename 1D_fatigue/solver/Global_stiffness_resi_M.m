function [KK_tan,F_gob]=Global_stiffness_resi_M(W,EL,d_n,phase_n1,infoPHF)
%%===================================
% calculates the stiffness matrix and internal force vector
%===================================
element=infoPHF.element;
eta=infoPHF.eta;
total_unknown=infoPHF.total_unknown;
numelem=infoPHF.numelem;
Ndof=infoPHF.Ndof; 
E=infoPHF.E;
Ce=E;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
K_vec = zeros(numelem*nnz_KE,1);
F_gob = zeros(total_unknown,1);
I_vec = zeros(numelem*nnz_KE,1);
J_vec = I_vec;
%===================================
for ele = 1 : numelem
    sctr=EL(ele,1).sctr;          % element connectivity
    nn = length(sctr) ;         % number of nodes per element
    disp_u=d_n(sctr,:);    % element disp at the respective nodes[u1 u2 u3 u4 u5 u6 u7 u8]
    phase=phase_n1(sctr); % element phase at the respective nodes [s1 s2 s3 s4]
    % -------------------------------------------------------
    % element stiffness and internal force vector
    K_e=zeros(nn*Ndof,nn*Ndof); f_int_e=zeros(nn*Ndof,1);
    % -------------------------------------------------------
    % Looping on Gauss point
    for kk = 1 : size(W,1)
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        N= EL(ele,1).GaussValues(kk,1).N;
        strain=Bfem*disp_u;
        %------------------------------------------------
        phase_gp=(N' * phase);
%       g_c=((1-eta)*(1-phase_gp)^(2))+eta;
        g_c=(1-phase_gp)^(2)+eta;  
        %------------------------------------------------
        [sigma,CC]= VM_2d_iso_M(strain,Ce,g_c);
        %-----------------------------------------------
        K_e=K_e+(Bfem'*CC*Bfem*W(kk)*detJ);
        f_int_e=f_int_e+(Bfem'*sigma*W(kk)*detJ);       
    end  % end of looping on GPs    
    I_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = repmat(sctr',num_DOF,1);
    tmp = repmat(sctr',1,num_DOF)';
    J_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = tmp(:);
    K_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = K_e(:);% Assemble the structural tangent element stiffness matrix and internal force vector
    F_gob(sctr,1)= F_gob(sctr,1) + f_int_e;
    
end % end of looping on elements

KK_tan = sparse(I_vec,J_vec,K_vec,total_unknown,total_unknown);

end
%=====================================================================

function[sigma,CC]= VM_2d_iso_M(strain,Ce,g_c)
%------------------------   
sigma = g_c*Ce*(strain);
CC=g_c*Ce;
%------------------------ 
end
%=====================================================================
