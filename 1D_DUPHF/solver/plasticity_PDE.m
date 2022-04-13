%================================================
% Plasticity
%================================================
% 2. solving the momentum equation
%------------------------------------------------
for l_it = 1 : niter
     [KK_plast_i,f_int_plast_i]=Global_plasticity_resi_M(W,EL,u_n1_i,...
            phase_n1_i,alpha_n_i,alpha_n,Eplas_n_i,infoPHF);
	%--------------------------------------------------
    C = KK_plast_i;
    dalpha = -C\f_int_plast_i; % solve for iterative increment of Det_u
    %--------------------------------------------------
	alpha_n_i = alpha_n_i + dalpha;       % update 
%------------------------------------------------
    % To determine plastic strain
    [Eplas_n_i]=cal_plastic_strain(W,EL,u_n1_i,phase_n1_i,...
        alpha_n_i,alpha_n,Eplas_n,infoPHF);    
    residual = norm(f_int_plast_i);      % residual 
   % disp(['PlastIt: ', num2str(l_it),' resi: ', num2str(residual)]);
    if (residual<tol)
        break
    end
end
%================================================
if l_it==niter && residual>0.1
   flag_terminate=1;
end
%------------------------------------------------
if infoPHF.irr_alpha == 1
	alpha_n_i = max(alpha_n_i,alpha_n);         % if ad hoc
end
%------------------------------------------------
% To determine plastic strain
[Eplas_n_i]=cal_plastic_strain(W,EL,u_n1_i,phase_n1_i,...
    alpha_n_i,alpha_n,Eplas_n,infoPHF);
%================================================
