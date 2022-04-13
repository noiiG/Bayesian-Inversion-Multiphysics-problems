function idof=node2dof(node,Ndof)

idof=[];

if Ndof==1
   for i=1:length(node)
       idof=[idof,node(i)];
   end 
elseif Ndof==2
   for i=1:length(node)
       idof=[idof,2*node(i)-1 2*node(i)];
   end
elseif Ndof==3
   for i=1:length(node)
       idof=[idof,3*node(i)-2 3*node(i)-1 3*node(i)];
   end
end   
