function [EL]=Global_resi_Initial(Q, W,elemType,infoPHF)

% calculates the stiffness matrix and internal force vector

%--------------------------------
element=infoPHF.element;
node=infoPHF.node;
total_unknown=infoPHF.total_unknown;
numelem=infoPHF.numelem;
element=infoPHF.element;
Ndof=infoPHF.Ndof; 
E=infoPHF.E;
%--------------------------------
num_DOF = size(element,2);
nnz_KE  = num_DOF^2;
%--------------------------------
for ele = 1 : numelem 
    sctr = element(ele,:);                % element connectivity
    nn = length(sctr) ;                   % number of nodes per element    
    EL(ele,1).sctr=sctr;
    % -------------------------------------------------------
    % B-operator (kinematics), element stiffness and internal force vector    
    % -------------------------------------------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                         % quadrature point
        % Shape functions and their derivatives
        [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx = dNdxi*invJ0;                        % derivatives of N, wrt XY
        B_phase=dNdx';   
        Bfem=B_phase ;% B matrix (corresponding to [u1 u2 u3 u4 u5 u6 u7 u8]) 
        %=================================================   
        % phase field terms     
        EL(ele,1).GaussValues(kk,1).N=N;
        EL(ele,1).GaussValues(kk,1).B=Bfem;
        EL(ele,1).GaussValues(kk,1).detJ0=det(J0);
        EL(ele,1).GaussValues(kk,1).B_phase=B_phase;
        Eplas(ele,1).data(kk,1).s=[0]';
    end  % end of looping on GPs
   %-----------------------------------------------------------      
end % end of looping on elements

