%================================================
% Initialaization for new step
%================================================
%disp('********************')
 %   disp([' LS: ', num2str(step)]);
    du_load_n = u_incr;
    du_n1(dof_disp) = du_load_n;                              % applied on the corresponding DOF
    u_n1(dof_disp) = u_n1(dof_disp)+du_load_n;                % increased U, not yet balanced (terms at the DOFs of Neumann BC)
    u_n1_i = u_n1;
    phase_n1_i=phase_n1;
    u_n1_0 = u_n1;    
    phase_n1_0=phase_n1;
%================================================