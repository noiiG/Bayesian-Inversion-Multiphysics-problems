namein=[namepre,'/'];
%====================
dl=lenB/nel;
node=linspace(0,lenB,nel+1)';
for i=1:nel
    element(i,:)=[i,i+1];
end
%nodeT=[node,zeros(1,length(node))];
%====================
[dof_nbc,dof_disp,total_unknown,numelem,numnode,dof_fix,dof_all_c]=BCGlobal(node,element,Ndof,lenB);