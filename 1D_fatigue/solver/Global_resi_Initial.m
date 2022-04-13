function [EL,gamma]=Global_resi_Initial(Q, W,elemType,infoPHF)

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
K_vec = zeros(numelem*nnz_KE,1);
%--------------------------------
for ele = 1 : numelem 
    sctr = element(ele,:);                % element connectivity
    nn = length(sctr) ;                   % number of nodes per element    
    EL(ele,1).sctr=sctr;
    % -------------------------------------------------------
    % B-operator (kinematics), element stiffness and internal force vector
    K_e=zeros(nn*Ndof,nn*Ndof);       
    % -------------------------------------------------------
    for kk = 1 : size(W,1)
        pt = Q(kk,:);                         % quadrature point
        % Shape functions and their derivatives
        [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx = dNdxi*invJ0;                        % derivatives of N, wrt XY
        B_phase=dNdx';
        % B matrix (corresponding to [u1 u2 u3 u4 u5 u6 u7 u8])        
        Bfem=B_phase ;
        tmp2 = Bfem'*E*W(kk)*det(J0);
        K_e=K_e+(tmp2*Bfem);
        %=================================================   
        % phase field terms     
        EL(ele,1).GaussValues(kk,1).N=N;
        EL(ele,1).GaussValues(kk,1).B=Bfem;
        EL(ele,1).GaussValues(kk,1).detJ0=det(J0);
        EL(ele,1).GaussValues(kk,1).B_phase=B_phase;
        gamma(ele,1).data(kk,1).s=[0]';
    end  % end of looping on GPs
    PI_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = repmat(sctr',num_DOF,1);
    tmp = repmat(sctr',1,num_DOF)';
    PJ_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = tmp(:);
    K_vec((ele-1)*nnz_KE+1:ele*nnz_KE) = K_e(:);
   %-----------------------------------------------------------      
end % end of looping on elements
KKn = sparse(PI_vec,PJ_vec,K_vec,total_unknown,total_unknown);

