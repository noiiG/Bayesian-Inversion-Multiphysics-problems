flag=0;
close all

figure
plot(results.MCMC(1,:),'LineWidth',1.5)
xlabel('iteration','FontWeight','norm','FontSize',55)
%ylabel('\alpha','FontWeight','norm','FontSize',55)
ylabel('G_c','FontWeight','norm','FontSize',55)

set(gca,'FontSize',50,'FontName','Times')

set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 8.5 6]);
%ylim([60000 74000]


 

figure
[y1,x1] = ksdensity(results.MCMC(1,1:1000),'Support','positive','Bandwidth',0.02);
plot(x1,y1,'LineWidth',3.5)
xlabel('G_c','FontWeight','norm','FontSize',55)
ylabel('pdf','FontWeight','norm','FontSize',55)
set(gca,'FontSize',50,'FontName','Times')
xline(0.04,'LineWidth',4,'LineStyle','--','Color','green');

%ylim([0 2.2e4])

set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 8.5 6]);
%ylim([0,16.5])


    

