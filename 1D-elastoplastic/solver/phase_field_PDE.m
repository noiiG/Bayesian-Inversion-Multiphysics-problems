%================================================
% Phase-field
%================================================
% 1. solving the phase field equation
%--------------------------------------------------
    set_cracked=find(phase_n1_i <= 0.01);
    dof_nbc_c = setdiff(dof_all_c, set_cracked);
%--------------------------------------------------
    [KK_phase,H_n1_new,F_phase,PI_vec,PJ_vec]=phase_field_iter_ISO_M_mex(EL, W,H_n1_i,u_n1_i,...
             Eplas_n_i,phase_n1_i,alpha_n_i,infoPHF);
            
    R_phase=F_phase(dof_nbc_c);          % residual r_phase
    %-----------------
    B = distributed(KK_phase(dof_nbc_c,dof_nbc_c));
    phase_n1_i(dof_nbc_c) = B\R_phase; % solve for iterative increment of Det_u
    %-----------------
    Min_c=min(phase_n1_i);
    if (Min_c<0)
    % disp('Warning: Neg phase values found at nodes !');
    Ind=find(phase_n1_i < 0);
    phase_n1_i(Ind)=0;
    end
%================================================    
