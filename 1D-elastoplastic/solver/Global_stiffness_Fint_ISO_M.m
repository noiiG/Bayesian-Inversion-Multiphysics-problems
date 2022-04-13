function [KK_tan,F_gob,alpha_n_i,Eplas_new]=Global_stiffness_Fint_ISO_M(W,EL,d_n,...
    alpha_n,Eplas,infoPHF)
%%===================================
alpha_n_i=alpha_n;
Eplas_new=Eplas;
%%===================================
% calculates the stiffness matrix and internal force vector
%===================================
total_unknown=infoPHF.total_unknown;
element=infoPHF.element;
numelem=infoPHF.numelem;
Ndof=infoPHF.Ndof; 
sigmaY=infoPHF.sigmaY;
E=infoPHF.E;
H=infoPHF.H;
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
    % -------------------------------------------------------
    % element stiffness and internal force vector
    K_e=zeros(nn*Ndof,nn*Ndof); f_int_e=zeros(nn*Ndof,1);
    % -------------------------------------------------------
    % Looping on Gauss point
    for kk = 1 : size(W,1)
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        N= EL(ele,1).GaussValues(kk,1).N;
        alpha0=alpha_n(ele,kk);
        strain_Pn=Eplas(ele,1).data(kk,1).s;
        strain=Bfem*disp_u;
        %------------------------------------------------
        [sigma,CC,alpha_i,strain_p]= VM_2d_iso_M(strain,strain_Pn,alpha0,H,sigmaY,E);
        alpha_n_i(ele,kk)=alpha_i;
        Eplas_new(ele,1).data(kk,1).s=strain_p;
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

function[sigma,CC,alpha,strain_p]= VM_2d_iso_M(strain,strain_Pn,alpha_i,H,sigmaY,E)
%------------------------
alphan=alpha_i;
qtrial=E*(strain-strain_Pn);
J2sigmaY = (sigmaY+(H*alphan));
PHI = abs(qtrial)-J2sigmaY;
%------------------------
if(PHI<=0)
    % ELASTIC STEP - UPDATE ELASTIC STRESS AND STRAINS
    alpha=alphan;
    strain_p=strain_Pn;
    sigma=qtrial;
    CC=E;
else
    % PLASTIC STEP    
    delgamma = (PHI)/(E+H); % LINEAR HARDENIN CASE
    % PLASTIC STEP   
    [normS]= cal_norm_J2S(strain,strain_Pn);
    alpha = alphan+(delgamma);  % UPDATE alpha 
    % Plastic STRAIN UPDATE
    strain_p=strain_Pn+delgamma*normS;
    % Elasto plastic tangent modulous
    sigma=qtrial-delgamma*E*normS;
    CC=E*H/(E+H);
end
end
%=========================================================================
function[snorm]= cal_norm_J2S(strain,strain_Pn)
%------------------------
snorm=sign(strain-strain_Pn);
end
%=========================================================================




























