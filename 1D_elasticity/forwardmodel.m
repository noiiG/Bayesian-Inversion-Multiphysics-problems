function [d1] = forwardmodel(e)
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]

%=========================================

%e=4.259;  % TO BE CHANGED
%=========================================
%Set material properties and length of the rod
a=1;
que=14;
etr=0.5;
l=1;
%=========================================
%FEM
nel=1000;
[d1,dl1,sig1,step1,nel1] = linelast(nel,e,a,que,etr,l);
%=========================================
 
 
end

