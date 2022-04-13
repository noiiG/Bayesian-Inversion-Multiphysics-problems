%================================================
% Elasticity
%================================================
% 2. solving the momentum equation
%------------------------------------------------
for l_it = 1 : niter
    [KK_n1_i,f_int_n1_i]=Global_stiffness_resi_M(W,EL,u_n1_i,...
            phase_n1_i,infoPHF);
    resi_nbc = -f_int_n1_i(dof_nbc);
    residual = norm(resi_nbc);
    %disp(['DispIt: ', num2str(l_it),' resi: ', num2str(residual)]);
    if (residual<tol)
        break
    end
    %-----------------------------------------------
    A = KK_n1_i(dof_nbc,dof_nbc);
    du_nbc = A\resi_nbc; % solve for iterative increment of Det_u
    %----------------------------------------------- increment of Det_u
    u_n1_i(dof_nbc) = u_n1_i(dof_nbc)+real(du_nbc);       % update
end
%================================================
if l_it==niter && residual>0.1
   flag_terminate=1;
end
%================================================
