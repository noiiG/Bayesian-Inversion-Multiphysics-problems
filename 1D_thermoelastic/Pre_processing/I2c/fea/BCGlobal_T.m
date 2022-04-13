%============================================================
% Boundary Conditions
%============================================================
function [dof_nbc,dof_disp,total_unknown,numelem,numnode,dof_fix,T_n]=BCGlobal_T(node,element,Ndof,lenB)

%========================================
%lenB=10;
%========================================
% Initilaization
%========================================
numnode = size(node,1); 
numelem = size(element,1);
total_unknown = numnode*Ndof;
dof_all = (1:1:total_unknown);
u_n = zeros(total_unknown,1); 
T_n=zeros(numnode,1);%30*ones(numnode,1); %zeros(numnode,1);
%========================================
% boundary conditions: Dirichlet BC and Neumann BC: DOF
%========================================

dof_fix=[];
dof_nbc=(2:1:numnode-1);
dof_disp=[1,numnode];
T_n(1)=1;
T_n(end)=0;





