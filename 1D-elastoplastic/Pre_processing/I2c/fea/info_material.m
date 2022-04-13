%=========================================
% Mtarix material properties
%=========================================
function info_matrirx=info_material(E,H,sigmaY,Ndof)
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
%=========================================
% Output Material Properties
%=========================================
info_matrirx.E=E;
info_matrirx.H=H;
info_matrirx.sigmaY=sigmaY;
info_matrirx.Ndof=Ndof;    
info_matrirx.C_vol=C_vol;
info_matrirx.C_dev=C_dev;
info_matrirx.Ce=Ce;
%================================================
