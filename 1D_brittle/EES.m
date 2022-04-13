function [ss1,d1] = EES(measurement,e)
 
[d1] = forwardmodel(e);

ss1=norm(measurement-d1);

end

