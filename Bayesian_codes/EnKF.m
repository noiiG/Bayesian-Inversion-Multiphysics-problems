function [results] = EnKF(inp)
% The EnKF-MCMC

format long

samples=inp.nsamples;
% start MCMC using MH algorithm tehn use ensemble Kalman
inp.nsamples=inp.Kalmans-1;
[results] = DRAM(inp);
disp('EnKF started')

% The obtained MCMC chain
X=results.MCMC;
Y=results.values;
thetaj=X(:,inp.Kalmans-1);
[oldpi,oldvalue]= EES(inp.measurement,thetaj);
accepted=fix(results.accepted*(inp.Kalmans-1)/100);

for j=inp.Kalmans:samples
    
    ss2=inp.me*ones(size(inp.measurement,1),1);
    RR = diag(ss2);
    
    
    mX = repmat(mean(X,2),1,j-1);
    mY = repmat(mean(Y,2),1,j-1);
    
    Ctm =  (X-mX)*(Y-mY)'/(j-2);
    Cmm =  (Y-mY)*(Y-mY)'/(j-2);
    
    % Kalman gain
    KK = Ctm/(Cmm+RR);
    
    [XX] = forwardmodel(thetaj);
    
    dt = KK*(inp.measurement+randn(length(inp.measurement),1).*inp.me-XX);
    thetas = thetaj+dt;

    
    thetas = check_bound(thetas,inp.range);
    
    [newpi,newvalue]=EES(inp.measurement,thetas);
    
    lambda = min(1,exp(-0.5*(newpi-oldpi)/inp.sigma ));
    
    % accept/reject the forecast
    if rand < lambda
        flag   = true;
        accepted     = accepted+1;
        thetaj   = thetas;
        oldpi    = newpi;
        oldvalue=newvalue;
    end
    
    X=[X,thetaj];
    Y=[Y,oldvalue];
%     lambda
     double(thetaj)
  
    
    if mod(j,100)==0
         fprintf('number of sample: %d\n',j)
         save('results','X','accepted')
     end
    
end

% The MCMC chain and the accetptance rate
results.MCMC=X;
results.accepted=(accepted/samples)*100;


end

