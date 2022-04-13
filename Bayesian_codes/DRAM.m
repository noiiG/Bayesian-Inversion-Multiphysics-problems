function [results] = DRAM(inp)
% DRAM algorithm

format long

epsilon=1e-5;

adaptation_number=100;

nsamples=inp.nsamples;
measurement=inp.measurement;

Rj = chol(inp.icov);

dim=size(inp.range,1);
Kp=2.4/sqrt(dim);


MCMC=zeros(nsamples,dim);
[oldpi,oldvalue]= EES(measurement,inp.theta0);
results.values=oldvalue;

 R2=Rj./5;
 invR=inv(Rj);

accepted=0;
MCMC(1,:)= inp.theta0;

thetaj=inp.theta0';


MCMCcov = []; MCMCmean = []; ss = []; ii = 0;

for j=2:nsamples
    
    flag = 0;
    thetas = thetaj+Rj*randn(dim,1);     % a new proposal
    
    
      thetas = check_bound(thetas,inp.range);
      
   
        [newpi,newvalue]=EES(measurement,thetas);
        lambda = min(1,exp(-0.5*(newpi-oldpi)/inp.sigma ));
        if rand < lambda % we accept
            flag   = true;
            accepted     = accepted+1;
            thetaj   = thetas;
            oldpi    = newpi;
            oldvalue=newvalue;
        end
        
        % lambda 
    
    if  ~flag 
        thetass = thetaj+R2*randn(dim,1);  % a new try
        
       thetass = check_bound(thetass,inp.range);
            
  
        [newss2,newvalue2]=EES(measurement,thetass);  
        k1 = min(1,exp(-0.5*(newpi-newss2)/inp.sigma ));
        k2 = exp(-0.5*(newss2-oldpi)/inp.sigma );
        k3 = exp(-0.5*(norm(invR*(thetass-thetas))^2));
         lambda2 = k2*k3*(1-k1)/(1-lambda);
         if rand < lambda2
             flag   = true;
             accepted     = accepted+1;
             thetaj   = thetass;
             oldpi    = newss2;
             oldvalue=newvalue2;
         end
      
    end
    
     results.values=[results.values,oldvalue];
    
     MCMC(j,:) = thetaj;
    % double(thetaj)
    
      if mod(j,adaptation_number)==0
         [MCMCcov,MCMCmean,ss] = covupd(MCMC((ii+1):j,:),1,MCMCcov,MCMCmean,ss);
         
         ii = j;
         [Ra,~] = chol(MCMCcov + eye(dim)*epsilon);

             Rj = Ra*Kp;
             
       
     end
    
    
      if mod(j,100)==0
         fprintf('number of sample: %d\n',j)
         save('results','MCMC','accepted')
     end
  
end

results.MCMC=MCMC';
results.accepted=(accepted/nsamples)*100;


end

