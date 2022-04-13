%================================================
% Phase-field
%================================================
% 1. solving the phase field equation
%--------------------------------------------------
set_cracked=find(phase_n1_i >= 0.999);
phase_n1_i(set_cracked)=1;
dof_nbc_c = setdiff(dof_all_c, set_cracked);
%--------------------------------------------------
for l_it = 1 : niter
    [KK_phase,H_n1_new,F_phase,PI_vec,PJ_vec]=phase_field_resi_M(EL, W,H_n1_i,phase_n1_i,u_n1_i,...
             phase_n1,ET_n_i,infoPHF); 	
    if condest(KK_phase(dof_nbc_c,dof_nbc_c)) <= 1e16    % quick check (note: inefficient)
		%--------------------------------------------------
        B = KK_phase(dof_nbc_c,dof_nbc_c);
        dphase = -B\F_phase(dof_nbc_c); % solve for iterative increment of Det_u
        %--------------------------------------------------
	else
		dphase = zeros(size(dof_nbc_c,1),1);
    end           
	phase_n1_i(dof_nbc_c)=phase_n1_i(dof_nbc_c) + dphase; 
    residual = norm(F_phase(dof_nbc_c));      % residual r_phase

  %  disp(['PhaseIt: ', num2str(l_it),' resi: ', num2str(residual)]);
    if (residual<tol)
        break
    end
  
end
%================================================   
if l_it==niter && residual>0.1
   flag_terminate=1;
end
%------------------------------------------------
if infoPHF.irr_phase == 1
	phase_n1_i = min(max(phase_n1_i,phase_n1),1);          % if ad hoc
end
%================================================    
