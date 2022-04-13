%============================================================
% Boundary Conditions
%============================================================
function [dof_nbc,dof_disp,total_unknown,numelem,numnode,dof_fix,dof_all_c]=BCGlobal(node,element,Ndof,lenB)

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
%========================================
% boundary conditions: Dirichlet BC and Neumann BC: DOF
%========================================
bound = find(node(:)==0)';            % left edge node fixed in   
force=find(node(:)==lenB)'; 

dof_fix0= bound;
dof_fix=[dof_fix0];
dof_disp= force;  

dof_nbc = setdiff(dof_all, [dof_fix,dof_disp]); % DOF need to solve/ active
%========================================
dof_all_c=(2:1:numnode-1);

