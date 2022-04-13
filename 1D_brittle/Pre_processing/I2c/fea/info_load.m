%================================================
% load info.
%================================================
% DISPLACEMENT controlled
u_incr = 1e-5;
nummulti=20;         % maximum number of iterations of
nstep=150;
niter1 = 5;                     % multiscale iteration  
% Newton-Raphson
niter = 30;           % number of Newton Raphson iterative steps
tol= 1e-8;          % tolerance of NR of local
no_conv=0;
niter_c=30;
no_conv_c=0;
opt_postmove=1;     % Post-processing
% Option for Composite Materials
optcomposite=0; % 1:=on 0:=off
% Option for non-matching mesh 
opt_nonmatch=1;  % 1:=on 0:=off

%================================================
% mesh type
%================================================
%--------------------
% Global
%--------------------
% Uniform meshing with Q4 elements
%elemType_G = 'T3' ;
%elemType = 'Q4' ;
elemType = 'L2' ;

% Gauss quadrature
%[W_G,Q_G] = quadrature_higher(1,'TRIANGULAR',1);
%[W,Q] = quadrature_higher(2,'GAUSS',2);
[W,Q] = quadrature_higher(1,'GAUSS',1);
%================================================

