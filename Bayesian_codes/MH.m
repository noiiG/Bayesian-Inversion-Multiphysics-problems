function [results] = MH(inp)
% The Metropolis-Hastings algorithm


nsamples=inp.nsamples;
measurement=inp.measurement;

Rj = chol(inp.icov);

% dimension of the parameters
dim=size(inp.range,1);

% set MCMC zero
MCMC=zeros(nsamples,dim);
[oldpi,oldvalue]= EES(measurement,inp.theta0);
results.values=oldvalue;

accepted=0;
MCMC(1,:)= inp.theta0;


thetaj=inp.theta0';


for j=2:nsamples
    
    flag = 0;
    % propose a new candidate
    thetas = thetaj+Rj*randn(dim,1);  
    
    
    % check the boundaries
    thetas = check_bound(thetas,inp.range);  
    
        [newpi,newvalue]=EES(measurement,thetas);
        lambda = min(1,exp(-0.5*(newpi-oldpi)/inp.sigma ));
       
        % accept/reject the candidate
        if rand < lambda  
            flag   = true;
            accepted     = accepted+1;
            thetaj   = thetas;
            oldpi    = newpi;
            oldvalue=newvalue;
        end
    
        results.values=[results.values,oldvalue];
    
        % update the MCMC chain
        MCMC(j,:) = thetaj;
      
        
      if mod(j,500)==0
         fprintf('number of sample: %d\n',j)
     end
  
end

results.MCMC=MCMC';
results.accepted=(accepted/nsamples)*100;


end

