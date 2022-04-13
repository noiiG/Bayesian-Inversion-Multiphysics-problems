% Initialization of: structural information
u_n1 = zeros(total_unknown,1);          % vector: nodal displacement
du_n1 = zeros(total_unknown,1);         % vector: increment of nodal disp.
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
inp_model.E=E;
inp_model.Ndof=Ndof;
inp_model.ls=ls;
inp_model.G_f=G_f;
inp_model.gamma_c=gamma_c;
inp_model.k=k;
%----------------------------------------
%phase_n1(1+(numnode-1)*0.5,1)=0.6;