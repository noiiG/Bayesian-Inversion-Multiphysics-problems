function [results] = MH_DR(inp)
% The delayed rejection Metropolis

nsamples=inp.nsamples;
measurement=inp.measurement;

Rj = chol(inp.icov);

dim=size(inp.range,1);

MCMC=zeros(nsamples,dim);
[oldpi,~]= EES(measurement,inp.theta0);

 R2=Rj./5;
 invR=inv(Rj);

accepted=0;
MCMC(1,:)= inp.theta0;

thetaj=inp.theta0';

for j=2:nsamples
    
    flag = 0;
    thetas = thetaj+Rj*randn(dim,1);     % a new proposal
    
    thetas = check_bound(thetas,inp.range);
        
   
        newpi=EES(measurement,thetas);
        lambda = min(1,exp(-0.5*(newpi-oldpi)/inp.sigma ));
        if rand < lambda % we accept
            flag   = true;
            accepted     = accepted+1;
            thetaj   = thetas;
            oldpi    = newpi;
        end
     
    
    if  ~flag 
        thetass = thetaj+R2*randn(dim,1);  % a new try
        
       thetass = check_bound(thetass,inp.range);
            
      
        [newss2,~]=EES(measurement,thetass);  
        k1 = min(1,exp(-0.5*(newpi-newss2)/inp.sigma ));
        k2 = exp(-0.5*(newss2-oldpi)/inp.sigma );
        k3 = exp(-0.5*(invR*norm((thetass-thetas))^2));
         lambda2 = k2*k3*(1-k1)/(1-lambda);
         if rand < lambda2
             flag   = true;
             accepted     = accepted+1;
             thetaj   = thetass;
             oldpi    = newss2;
        
        
        end
    end
    
    
    
    MCMC(j,:) = thetaj;
    
      if mod(j,500)==0
         fprintf('number of sample: %d\n',j)
     end
  
end

results.MCMC=MCMC';
results.accepted=(accepted/nsamples)*100;


end

