%=========================================
% Private FEM code for 2D linear problems
%=========================================
% @copyright(2020)
% This code is for ELASTICITY written by
% ** Nima Noii **
% Institute of Contiunum Mechancis 
% (IKM-LUH, Hannover, DE)
%=========================================
%
function [d,dl,sig,step,nel] = linelast(nel,e,a,que,etr,l)
ea=e*a;
%Set up number of equations
ng=nel+1;

%Initialize global displacement vector for initial loop
d=zeros(nel+1,1);

%Compute length of one element
dl=l/nel;

%SET UP ELEMENT MATRICES
%Element stiffness matrix
ke=ea/dl*[1 -1;-1 1];

%SET UP GLOBAL MATRICES
%Initialize load vector and global stiffness matrix to zero
q=zeros(ng,1);
kg=zeros(ng);

%Assemble the global stiffness and load vector
for i=1:nel
  x1=(i-1)*dl;
  x2=i*dl;
  for j=1:2
    for k=1:2
      kg(j+i-1,k+i-1)=kg(j+i-1,k+i-1)+ke(j,k);
    end
    if j==1
    q(i)=q(i)+(l*cos(2.0*pi*x1/l)/(2.0*pi)... 
             -l*l*sin(2.0*pi*x2/l)/(4.0*pi*pi*dl)...
	       +l*l*sin(2.0*pi*x1/l)/(4.0*pi*pi*dl))*que;
    elseif j==2
    q(i+1)=q(i+1)+(-l*cos(2.0*pi*x2/l)/(2.0*pi)... 
                 +l*l*sin(2.0*pi*x2/l)/(4.0*pi*pi*dl)...
		     -l*l*sin(2.0*pi*x1/l)/(4.0*pi*pi*dl))*que;
    end
  end
end

%BOUNDARY CONDITIONS AND SOLUTION
%Put essential boundary conditions d(x=0)=0
kg(1,:)=[1 zeros(1,nel)];
kg(:,1)=[1; zeros(nel,1)];
q(1)=0;

%Add applied traction to the last node
q(ng)=q(ng)+etr*a;

%Compute displacement
d(1:ng,1)=kg\q;  % kg^-1*q

%Compute stresses: post-processing
for i=1:nel
  sig(2*i-1)=e/dl*(d(i+1)-d(i));
  sig(2*i)=e/dl*(d(i+1)-d(i));
  step(2*i-1)=dl*(i-1);
  step(2*i)=dl*i;  
end
