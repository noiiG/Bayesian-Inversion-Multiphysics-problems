%=========================================
% Mtarix material properties
%=========================================
function info_matrirx=info_material(E,Ndof,G_f,irr_phase,ls,Calpha,rho_C_p,Ck)
%=========================================
% Mtarix material properties
%=========================================
%=========================================
%----------------------------------------
% Elasticity
%----------------------------------------
Ce = E;
C_vol =[];
C_dev =[];
%----------------------------------------
% Phase-field
%----------------------------------------
%ls=0.3;    % =2*ls 
eta=1e-6;    % small dimensionless parameter 0< eta <<1
Fact1=4*(ls^2);
Fact2=4*ls*(1-eta);
%=========================================
% Output Material Properties
%=========================================
info_matrirx.E=E;
info_matrirx.G_f=G_f;
info_matrirx.irr_phase=irr_phase;

info_matrirx.eta=eta;
info_matrirx.Ndof=Ndof;    
info_matrirx.ls=ls;

info_matrirx.Fact1=Fact1;
info_matrirx.Fact2=Fact2;

info_matrirx.Calpha=Calpha;
info_matrirx.rho_C_p=rho_C_p;
info_matrirx.Ck=Ck;
%================================================
