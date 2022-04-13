function out_inter_refine=initial_info_inter

%=======================================================                            
% Initial Global Element that has to be refined 
%=======================================================  
IDeleR_refine=[39:90];         % Reference 
IDeleG_refine=[39:90];  % MultiScale

%--------------
% Output
%--------------
out_inter_refine.IDeleR_refine=IDeleR_refine;
out_inter_refine.IDeleG_refine=IDeleG_refine;