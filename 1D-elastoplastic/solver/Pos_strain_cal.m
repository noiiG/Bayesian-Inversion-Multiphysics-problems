function [Strain_pos,EP_bar]=Pos_strain_cal(strain_p)


strain_Pmat = [strain_p(1,1)    strain_p(3,1)/2 0;
               strain_p(3,1)/2  strain_p(2,1)   0;
               0                0               strain_p(4,1)];    % tensor form
                   
[Vec, Eig_v]=eig(strain_Pmat);
p_strain1=Eig_v(1,1);p_strain2=Eig_v(2,2);p_strain3=Eig_v(3,3);
    
p_strain1_P=(p_strain1+abs(p_strain1))/2;
p_strain2_P=(p_strain2+abs(p_strain2))/2;
p_strain3_P=(p_strain3+abs(p_strain3))/2;

Pos_eig=[p_strain1_P 0           0;
         0           p_strain2_P 0;
         0           0            p_strain3_P];

Strain_pos=Vec*Pos_eig*Vec';

E_P2=Strain_pos*Strain_pos;
tr_E_P=E_P2(1,1)+E_P2(2,2)+E_P2(3,3);
EP_bar=sqrt(2/3)*sqrt(tr_E_P);