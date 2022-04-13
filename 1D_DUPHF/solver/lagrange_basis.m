function [Nv,dNdxi]=lagrange_basis(type,coord,dim)
  
% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' 
%  
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3  
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University    
    
  if ( nargin == 2 )
    dim=1;
  end
  
  switch type
   case 'L2'  
    %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2
    %
    if size(coord,2) < 1
      disp('Error coordinate needed for the L2 element')
    else
     xi=coord(1);
     N=([1-xi,1+xi]/2)';
     dNdxi=[-1;1]/2;
    end
  
   case 'L3' 
    %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2----------3
    %
    if size(coord,2) < 1
      disp('Error two coordinates needed for the L3 element')
    else
     xi=coord(1);
     N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2];
     dNdxi=[xi-.5;xi+.5;-2*xi];
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
  end
 
  I=eye(dim);
  Nv=[];
  for i=1:size(N,1)
    Nv=[Nv;I*N(i)];
  end
  
  if ( dim == 1 )
    B=dNdxi;
    
  end  
  
  % end of function
