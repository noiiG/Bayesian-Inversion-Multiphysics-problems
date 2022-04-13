%================================================
% Elasticity
%================================================
% 2. solving the momentum equation
%------------------------------------------------
for l_it = 1 : niter
    [KK_n1_i,f_int_n1_i,alpha_n_i,Eplas_n_i]=Global_stiffness_Fint_ISO_M(W,EL,u_n1_i,...
            alpha_n,Eplas_n,infoPHF);
    resi_nbc = -f_int_n1_i(dof_nbc);
    residual = norm(resi_nbc);
   % disp(['It: ', num2str(l_it),' resi: ', num2str(residual)]);
    if (residual<tol)
        break
    end
    %-----------------------------------------------
    A = KK_n1_i(dof_nbc,dof_nbc);
    du_nbc = A\resi_nbc; % solve for iterative increment of Det_u
    %-----------------------------------------------
    u_n1_i(dof_nbc) = u_n1_i(dof_nbc)+real(du_nbc);       % update
    if l_it==niter && residual>0.1
       flag_terminate=1;
    end
end
alpha_n_new=alpha_n_i;
%================================================
