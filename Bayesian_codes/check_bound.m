function x = check_bound(x,range)
% check the boundary of the proposal according to the given bounds
% if it exceeds, set min or max


mini=range(:,1);  maxi=range(:,2);

imi=find(x<mini);   ima=find(x>maxi);

if ~isempty(imi)
x(imi)=mini(imi);
end

if ~isempty(ima) 
x(ima)=maxi(ima);
end 

end