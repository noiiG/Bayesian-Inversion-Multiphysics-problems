%================================================
% Initialaization for new step
%================================================
%disp('********************')
 %   disp([' LS: ', num2str(step)]);
    du_load_n = u_incr;
    du_n1(dof_disp) = du_load_n;                              % applied on the corresponding DOF
    u_n1(dof_disp) = u_n1(dof_disp)+du_load_n;                % increased U, not yet balanced (terms at the DOFs of Neumann BC)
    u_n1_i = u_n1;
    Eplas_n_i=Eplas_n;
    alpha_n_new=alpha_n;
     %----------------------
     % Preliminary calculation for the new U, using KK_n
     % disp('    Preliminary calculation');
     df_int_n1_pre = KK_n1_i*du_n1;
     resi_nbc = -df_int_n1_pre(dof_nbc);                      % residual r_nbc
     du_nbc = KK_n1_i(dof_nbc,dof_nbc)\resi_nbc;                 % solve for increment: du_nbc (since du_dbc=0)
     du_n1_tr = zeros(total_unknown,1);
     du_n1_tr(dof_nbc) = du_nbc;
     u_n1_i = u_n1_i+du_n1_tr;
     %----------------------
%================================================