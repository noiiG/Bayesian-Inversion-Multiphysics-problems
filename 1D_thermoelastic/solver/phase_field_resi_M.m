function [KK_phase,H_n1_new,F_phase,PI_vec,PJ_vec]=phase_field_resi_M(EL, W,...
    H_se_n,phase_n1_i,u_n1_i,phase_n1,ET_n_i,infoPHF)

% calculates the stiffness matrix and internal force vector

%%===================================
element=infoPHF.element;
G_f=infoPHF.G_f;
eta=infoPHF.eta;
E=infoPHF.E;
numelem=infoPHF.numelem;
numnode=infoPHF.numnode;
ls=infoPHF.ls;
node=infoPHF.node;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
PI_vec = zeros(numelem*nnz_KE,1);
PJ_vec = PI_vec;
K_vec = zeros(numelem*nnz_KE,1);
H_n1_new = zeros(numelem, length(W)); %  History variable
F_phase = zeros(numnode,1);
D=0;
drive=0;
drivejac=0;
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele).sctr;          % element connectivity
    nn = length(sctr) ;         % number of nodes per element
    disp_u=u_n1_i(sctr,:);     %element disp at the respective nodes[u1 u2 u3 u4 u5 u6 u7 u8]
    phase=phase_n1_i(sctr);       % current phf estimate 
    phase_n=phase_n1(sctr);       % current phf estimate 
    % -------------------------------------------------------
    K_p_e=zeros(nn,nn);
    f_int_p=zeros(length(EL(ele,1).sctr'),1);
    % ---------------------------------------------
    [A_aniso]=1;
    % ---------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        B_phase=EL(ele,1).GaussValues(kk,1).B_phase;
        % calculating field values at GP at current iteration
        strain=Bfem*disp_u;
        phase_gp=(N' * phase);
        phasen_gp=(N' * phase_n);
        Grad_phase=B_phase*phase;
%         diff_g=2*(1-eta)*(1-phase_gp);  
%         diff2_g=-2*(1-eta);
        diff_g=2*(1-phase_gp);
        diff2_g=-2;
        Fact1=ls^2;    
        %--------------------------------------
        % calculating driving force
        strain_T=ET_n_i(ele,1).data(kk,1).s;
        strain_e=strain-strain_T;
        strain_energy=strain_energy_ductile_iso_M(strain_e,E);    
        %--------------------------------------
        % calculating history variable
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
        % irreversibility
		D=ls*H_se/G_f;
		drive=diff_g*D-phase_gp;
		drivejac=diff2_g*D-1;
        %--------------------------------------
        % Jacobian
		K_phase = (N*(drivejac)*N')-(B_phase' *Fact1*A_aniso* B_phase);
		%--------------------------------------
		% Residual
		f_micro_force1= N * (drive)-(B_phase' *Fact1*A_aniso* Grad_phase);
		K_p_e=K_p_e+(K_phase*W(kk)*detJ);
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
 function strain_energy=strain_energy_ductile_iso_M(strain,E)      

%-------------------------------
strain_energy_ele=0.5*E*strain^2;
%-------------------------------
strain_energy=strain_energy_ele;
end
%=========================================================================



































