function gamma=cal_fat_field(W,EL,u_n1_i,u_n1,phase_n1_i,phase_n1,gamma,infoPHF)

%===================================
E=infoPHF.E;
numelem=infoPHF.numelem;
%===================================
for ele = 1 : numelem
    % ---------------------------------------------
    sctr=EL(ele,1).sctr;      % element connectivity
    
    disp_u=u_n1_i(sctr,:);    % element disp at the respective nodes
    disp_u_n=u_n1(sctr,:);    

    phase=phase_n1_i(sctr);   % element phase field at the respective nodes 
    phase_n=phase_n1(sctr); 
    
    % -------------------------------------------------------
    for kk = 1 : size(W,1)      % Looping on Gauss point
        Bfem= EL(ele,1).GaussValues(kk,1).B;
        %--------------------------------------
        % calculating field values at the GP 
        strain=Bfem*disp_u;
        strain_n=Bfem*disp_u_n;        
        %--------------------------------------
        phase=(phase(1)+phase(2))/2;
        phase_n=(phase_n(1)+phase_n(2))/2;
        psi_e=(1-phase)^2*0.5*E*(strain)^2;
        psi_e_n=(1-phase_n)^2*0.5*E*(strain_n)^2;
        
        if (psi_e-psi_e_n)>0
            gamma(ele,1).data(kk,1).s=gamma(ele,1).data(kk,1).s+(psi_e-psi_e_n);    
        end
        %--------------------------------------
    end  % end of looping on GPs
end % end of looping on elements
end






















