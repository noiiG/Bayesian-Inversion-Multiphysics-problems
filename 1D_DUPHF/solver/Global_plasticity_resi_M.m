function [KK_phase,F_phase]=Global_plasticity_resi_M(W,EL,d_n,phase_n1,...
    alpha_n_i,alpha_n,Eplas,infoPHF)

%%===================================
E=infoPHF.E;
numelem=infoPHF.numelem;
element=infoPHF.element;
numnode=infoPHF.numnode;
H=infoPHF.H;
sigmaY=infoPHF.sigmaY;
l_p=infoPHF.l_p;
irr_alpha=infoPHF.irr_alpha;
eta=infoPHF.eta;
%===================================
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
PI_vec = zeros(numelem*nnz_KE,1);
PJ_vec = PI_vec;
K_vec = zeros(numelem*nnz_KE,1);
F_phase = zeros(numnode,1);
penalty=0;
penaltygrad=0;
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele,1).sctr;          % element connectivity
    nn = length(sctr) ;         % number of nodes per element
    disp_u=d_n(sctr,:);    % element disp at the respective nodes[u1 u2 u3 u4 u5 u6 u7 u8]
    phase=phase_n1(sctr); % element phase at the respective nodes [s1 s2 s3 s4]
	alphap=alpha_n_i(sctr);
    alphapn=alpha_n(sctr);
    % -------------------------------------------------------
    K_p_e=zeros(nn,nn);
    f_int_p=zeros(length(EL(ele,1).sctr'),1);
    % ---------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        N= EL(ele,1).GaussValues(kk,1).N;
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        detJ=EL(ele,1).GaussValues(kk,1).detJ0;
        B_phase=EL(ele,1).GaussValues(kk,1).B_phase;
        strain_Pn=Eplas(ele,1).data(kk,1).s;
        % calculating field values at GP at current iteration
        strain=Bfem*disp_u;
        %--------------------------------------
        phase_gp=(N' * phase);
        alphapn_gp=(N' * alphapn);
		alphap_gp=(N' * alphap);
		Grad_alpha=B_phase*alphap;
        %--------------------------------------
        % irreversibility
		switch irr_alpha
		case 1
			penalty=0;
			penaltygrad=0;
		case 2
			pen_alpha_irr=1e9;           % irreversibility penalty
            penalty=-pen_alpha_irr*min(alphap_gp-alphapn_gp,0.0);
            if ((alphap_gp-alphapn_gp)<=0.0) 
				penaltygrad=-pen_alpha_irr;
			else			
				penaltygrad=0;
            end
        end
        %--------------------------------------
        g_c=(1-phase_gp)^(2)+eta;
        E_d=g_c*E;
        Y_d=g_c*sigmaY;
        H_d=g_c*H;
        %--------------------------------------
        % internal force: phase field
        [normS]= cal_J2S(strain,strain_Pn);
        strain_e=strain-strain_Pn;
        f_micro_force1=N*(-E_d*strain_e+H_d*alphap_gp+Y_d)+B_phase'*(Y_d*l_p^2)* Grad_alpha;
        f_int_p=f_int_p+ (f_micro_force1*W(kk)*detJ);
        %--------------------------------------    
        K_phase=N*(normS*E_d+H_d)*N'+B_phase'*(Y_d*l_p^2)* B_phase;
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
function[snorm]= cal_J2S(strain,strain_Pn)
%------------------------
snorm=sign(strain-strain_Pn);
end
%=========================================================================









