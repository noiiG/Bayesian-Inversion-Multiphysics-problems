%=========================================
function [Plot_P]=forwardmodel(x)
%=========================================
% Material and model parameters

% [ 1] [ E       ] Young's modulus ........................[MPa]
% [ 2] [ G_c     ] fracture toughness .....................[MPa mm]
% [ 3] [ l_f     ] fracture length scale ..................[mm]
% [ 4] [ gamma_c ] fatigue threshold ......................[MPa]
% [ 5] [ k       ] fatigue degradation parameter ..........[-]

%=========================================
psi_c=30;
Ndof=1;
 
%-----------------
lenB=1;
nel=300;
h=lenB/nel;
%-----------------
 

E=x(1);
k=x(2);
G_f=x(3);
gamma_c=x(4);
ls=0.55/(sqrt(2)*sqrt(psi_c));%2*h;   


%=========================================
[~,Plot_P]=main(E,Ndof,G_f,ls,lenB,nel,gamma_c,k);
%=========================================

end

