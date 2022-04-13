%================================================
% Plasticity
%================================================
% 2. solving the momentum equation
%------------------------------------------------
for l_it = 1 : niter
     [KK_thermo_i,f_int_thermo_i]=Global_stiffness_thermal_resi_M(W,EL,phase_n1_i,T_n_i,T_n,infoPHF);
    %--------------------------------------------------
    C = KK_thermo_i(dof_nbc_T,dof_nbc_T);
    dthet = -C\f_int_thermo_i(dof_nbc_T); % solve for iterative increment of Det_u
    %--------------------------------------------------
	T_n_i(dof_nbc_T) = T_n_i(dof_nbc_T) + dthet;       % update 
%------------------------------------------------   
    residual = norm(f_int_thermo_i(dof_nbc_T));      % residual 
    %disp(['ThermoIt: ', num2str(l_it),' resi: ', num2str(residual)]);
    if (residual<tol)
        break
    end
end
%================================================
if l_it==niter && residual>0.1
   flag_terminate=1;
end
%------------------------------------------------
% To determine plastic strain
[ET_n_i]=cal_thet_strain(W,EL,T_n_i,infoPHF);
%================================================
