function X = prior(N,Nx,range)
% The prior values

xmin = range(:,1);
xmax = range(:,2);

X = nan(N,Nx);

for i = 1:N
    X(i,:) = unifrnd(xmin,xmax);
end    

end