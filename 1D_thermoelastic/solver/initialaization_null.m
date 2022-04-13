% Initialization of: structural information
u_n1 = zeros(total_unknown,1);          % vector: nodal displacement
du_n1 = zeros(total_unknown,1);         % vector: increment of nodal disp.
T_n=T_n0;
T_n_i=T_n;
% Initialization of: phase field information
P_ext_n1=zeros(numnode,1);
phase_n1 = zeros(numnode,1);
H_n1_i=zeros(numelem,size(W,1));
p_n=zeros(numelem,size(W,1));
%--------------------------------
% initial calculations
I_vec  = []; J_vec = [];PI_vec  = []; PJ_vec = [];
%--------------------------------
Plot_P = zeros(nstep+1,1);
Plot_U = zeros(nstep+1,1);
flag_terminate=0;
%----------------------------------------
%phase_n1(1+(numnode-1)*0.5,1)=0.6;