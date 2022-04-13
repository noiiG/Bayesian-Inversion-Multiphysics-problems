function [W,Q] = quadrature_higher( quadorder, qt, sdim )

% The function quadrature returns a n x 1 column vector W of quadrature
% weights and a n x dim matrix of quadrature points, where n is the
% number of quadrature points.  The function is called as follows:
%
%   [W,Q]=quadrature( nint, type, dim )
%
% nint is the quadrature order, type is the type of quadrature
% (i.e. gaussian, triangular, etc.. ) and dim is the number of spacial
% dimentions of the problem.  The default for type is GAUSS and the
% default for dim is unity.
%
% wrQ=quadrature(nint,'TRIANGULAR',2);itten by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University


  if ( nargin < 3 )   % set default arguments
    if ( strcmp(qt,'GAUSS') == 1 )
      dim = 1;
    else
      dim = 2;
    end
  end

  if ( nargin < 2 )
    type = 'GAUSS';
  end

  if ( strcmp(qt,'GAUSS') == 1 ) 

    if ( quadorder > 8 )  % check for valid quadrature order
      disp('Order of quadrature too high for Gaussian Quadrature'); 
      quadorder =8;
    end
    
    quadpoint=zeros(quadorder^sdim ,sdim);
    quadweight=zeros(quadorder^sdim,1);
  
    r1pt=zeros(quadorder,1); r1wt=zeros(quadorder,1);

    switch ( quadorder ) 
      case 1
        r1pt(1) = 0.000000000000000;
        r1wt(1) = 2.000000000000000;

      case 2
        r1pt(1) = 0.577350269189626;
        r1pt(2) =-0.577350269189626;

        r1wt(1) = 1.000000000000000; 
        r1wt(2) = 1.000000000000000;         

      otherwise
        disp('Order of quadrature to high for Gaussian Quadrature'); 
	
    end  % end of quadorder switch

    n=1;
     
    if ( sdim == 1 ) 
      for i = 1:quadorder
        quadpoint(n,:) = [ r1pt(i) ];           
        quadweight(n) = r1wt(i); 
        n = n+1;
      end
      
    end
    
    Q=quadpoint;
    W=quadweight;
  % END OF GAUSSIAN QUADRATURE DEFINITION
  
    
  end  % end of TRIANGULAR initialization
  
% END OF FUNCTION
