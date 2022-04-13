function [KK_phase,H_n1_new,F_phase,PI_vec,PJ_vec]=phase_field_iter_ISO_M(EL, W,...
    H_se_n,u_n1_i,Eplas_n_i,phase_n_i,alpha_n,infoPHF)

% calculates the stiffness matrix and internal force vector
coder.allowpcode('plain')
%===================================
eta=infoPHF.eta;
lamda=infoPHF.lamda;
mu=infoPHF.mu;
numelem=infoPHF.numelem;
numnode=infoPHF.numnode;
EP_bar_crit=infoPHF.EP_bar_crit;
Fact1=infoPHF.Fact1;
Fact2=infoPHF.Fact2;
D_with_thresh=1;
psi_c=infoPHF.psi_c;
H=infoPHF.H;
sigmaY=infoPHF.sigmaY;
ls=infoPHF.ls;
chi=infoPHF.chi;
node=infoPHF.node;
theta=infoPHF.theta;
%===================================
num_DOF = 4;
nnz_KE  = num_DOF^2;

num_DOF = 4; nnz_KE  = num_DOF^2;
PI_vec = zeros(numelem*nnz_KE,1);
PJ_vec = PI_vec;
K_vec = zeros(numelem*nnz_KE,1);
H_n1_new = zeros(numelem, 4); %  History variable
F_phase = zeros(numnode,1);

for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele).sctr;          % element connectivity
    sctrB=EL(ele).sctrB;        % Scatter vector for element assembly (elemental DOFs)
    nn = length(sctr) ;         % number of nodes per element
    disp_u=u_n1_i(sctrB,:);    %element disp at the respective nodes[u1 u2 u3 u4 u5 u6 u7 u8]
    phase=phase_n_i(sctr);
    K_p_e=zeros(nn,nn);
    f_int_p=zeros(length(EL(ele,1).sctr'),1);
    % ---------------------------------------------
    [A_aniso]=fun_info_anisoM(chi,theta);
    % ---------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        B_phase=EL(ele,1).GaussValues(kk,1).B_phase;
        % calculating field values at GP at current iteration
        strain=Bfem*disp_u;
        strain_p=Eplas_n_i(ele,1).data(kk,1).s;
        strain_e=[strain-strain_p(1:3,1);-strain_p(4,1)];
        alpha=alpha_n(ele,kk);
        %--------------------------------------
        % calculating history variable
        strain_energy=strain_energy_ductile_iso_M(strain_e,H,alpha,sigmaY,lamda,mu);    
        if(strain_energy-H_se_n(ele,kk) > 0)
            H_se=strain_energy;
        else
            H_se=H_se_n(ele,kk);
        end
        if (H_se<0)
            disp('Warning: Neg H_se found !');
        end
        H_n1_new(ele,kk) = H_se;
        %--------------------------------------
        % phase field term
        Fact1=ls^2;
        if D_with_thresh==0     % Without Threshold
           D=ls*H_se/G_f;
        elseif D_with_thresh==1 % With Threshold
           D=0.5*fun_heavyside(H_se/psi_c-1)*(H_se/psi_c-1);
        end
        diff_g=2*(1-eta);      
        L_Fact=diff_g*D;
        K_phase = (N*((L_Fact)+1)* N')+(B_phase' *Fact1*A_aniso* B_phase);
        K_p_e=K_p_e+(K_phase*W(kk)*detJ);
        %--------------------------------------
        % internal force: phase field
        f_micro_force1= N * 1;
        f_int_p=f_int_p+ (f_micro_force1*W(kk)*detJ);
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
function Hout=fun_heavyside(Hin)
if Hin>0
    Hout=1;
else
    Hout=0;
end
end
%=========================================================================
function [A_aniso]=fun_info_anisoM(chi,theta)

% theta=30;

a_aniso=[cosd(theta);sind(theta)];
M_aniso=a_aniso*a_aniso';
A_aniso=eye(2)+chi*M_aniso;

end
%=========================================================================
 function strain_energy=strain_energy_ductile_iso_M(strain,H,alpha,sigmaY,lamda,mu)      


strain_mat = [strain(1) strain(3)/2 0;
                      strain(3)/2 strain(2) 0;
                       0         0           strain(4)];    % tensor form

E_trace = (strain_mat(1,1)+strain_mat(2,2)+strain_mat(3,3));

E_dev2=strain_mat*strain_mat;tr_E_dev=E_dev2(1,1)+E_dev2(2,2)+E_dev2(3,3);

strain_energy_ele=(E_trace^2*(lamda/2))+(mu*tr_E_dev);
strain_energy_pl=alpha*(sigmaY+0.5*H*alpha);
%-------------------------------
strain_energy=strain_energy_ele+strain_energy_pl;
end
%=========================================================================

