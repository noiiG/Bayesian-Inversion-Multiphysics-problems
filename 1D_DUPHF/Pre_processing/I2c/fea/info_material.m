%=========================================
% Mtarix material properties
%=========================================
function info_matrirx=info_material(E,psi_c,H,sigmaY,EP_bar_crit,p_stable,Ndof,chi,l_p,G_f,at,irr_phase,irr_alpha,ls,theta)
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
info_matrirx.at=at;
info_matrirx.irr_phase=irr_phase;
info_matrirx.irr_alpha=irr_alpha;
info_matrirx.H=H;
info_matrirx.sigmaY=sigmaY;
info_matrirx.EP_bar_crit=EP_bar_crit;

info_matrirx.eta=eta;
info_matrirx.Ndof=Ndof;    
info_matrirx.ls=ls;

info_matrirx.C_vol=C_vol;
info_matrirx.C_dev=C_dev;
info_matrirx.Ce=Ce;
info_matrirx.Fact1=Fact1;
info_matrirx.Fact2=Fact2;
info_matrirx.psi_c=psi_c;
%================================================
