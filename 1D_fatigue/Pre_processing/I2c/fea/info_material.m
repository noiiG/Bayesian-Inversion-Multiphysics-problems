%=========================================
% Mtarix material properties
%=========================================
function info_matrirx=info_material(E,Ndof,G_f,ls,gamma_c,k)
%=========================================
% Mtarix material properties
%=========================================
%=========================================
%----------------------------------------
% Elasticity
%----------------------------------------
Ce = [];
C_vol =[];
C_dev =[];
%----------------------------------------
% Phase-field
%----------------------------------------
%ls=0.3;    % =2*ls 
eta=1e-6;    % small dimensionless parameter 0< eta <<1
% eta=0;    % small dimensionless parameter 0< eta <<1
Fact1=4*(ls^2);
Fact2=4*ls*(1-eta);
%=========================================
% Output Material Properties
%=========================================
info_matrirx.E=E;
info_matrirx.G_f=G_f;
info_matrirx.gamma_c=gamma_c;
info_matrirx.k=k;


info_matrirx.eta=eta;
info_matrirx.Ndof=Ndof;    
info_matrirx.ls=ls;

info_matrirx.C_vol=C_vol;
info_matrirx.C_dev=C_dev;
info_matrirx.Ce=Ce;
info_matrirx.Fact1=Fact1;
info_matrirx.Fact2=Fact2;
%================================================
