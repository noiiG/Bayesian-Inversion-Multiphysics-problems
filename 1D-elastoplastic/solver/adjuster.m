function [out] = adjuster(inp,indi)

kk=4;
nn=size(inp,2);

for i=1:length(indi)

     tt=indi(i);
    if (tt+kk)>nn
        break
    end
    
   
    r1=inp(tt:tt+kk);
    
    b1=max(r1);
    a1=min(r1);
    
    r = ((b1-a1).*rand(kk,1) + a1)';
    
    inp=[inp(1:tt) r inp(tt+1:end)];
   
%     if i~=length(indi)
%     indi(i+1)=indi(i+1)+tt;
%     end
    
    %inp=inp(1:nn);
    
end

out=inp;

end

