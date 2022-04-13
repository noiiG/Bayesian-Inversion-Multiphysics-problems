function [results] = AMH(inp)
% adaptive Mtropois-Hastings algorithm


epsilon=1e-5;

adaptation_number=100;

nsamples=inp.nsamples;
measurement=inp.measurement;

Rj = chol(inp.icov);

dim=size(inp.range,1);
Kp=2.4/sqrt(dim);


MCMC=zeros(nsamples,dim);
[oldpi,~]= EES(measurement,inp.theta0);

accepted=0;
MCMC(1,:)= inp.theta0;

thetaj=inp.theta0';

MCMCcov = []; MCMCmean = []; ss = []; ii = 0;

for j=2:nsamples
    
    flag = 0;
    thetas = thetaj+Rj*randn(dim,1);     % a new proposal
    
    thetas = check_bound(thetas,inp.range); 
        
   
    [newpi,~]=EES(measurement,thetas);
    lambda = min(1,exp(-0.5*(newpi-oldpi)/inp.sigma ));
        if rand < lambda % we accept
            flag   = true;
            accepted     = accepted+1;
            thetaj   = thetas;
            oldpi    = newpi;
        end
   
    
    MCMC(j,:) = thetaj;
    
    if mod(j,adaptation_number)==0
         [MCMCcov,MCMCmean,ss] = covupd(MCMC((ii+1):j,:),1,MCMCcov,MCMCmean,ss);
         
         ii = j;
         [Ra,~] = chol(MCMCcov + eye(dim)*epsilon);

             Rj = Ra*Kp;
             
       
     end
    
      if mod(j,500)==0
         fprintf('number of sample: %d\n',j)
     end
  
end

results.MCMC=MCMC';
results.accepted=(accepted/nsamples)*100;


end

