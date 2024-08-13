function plot_ev(betavalue,timelab,titlename,ylab,setcolor,ifplot,nplot)
% output: a plot with ev beta means, std error and significant points marked
% input: betavalue,timelab,titlename,setcolor
%
if nargin<6
    ifplot=true;
end
if nargin<7
    nplot=1;
end
mean1=squeeze(mean(betavalue))';
sem=squeeze(std(betavalue)./sqrt(size(betavalue,2)));
sem=sem';
jbfill(timelab,mean1+sem,mean1-sem,setcolor,setcolor,0,0.2);

hold on
plot(timelab,mean1,'LineWidth',2,'color',setcolor);
hold on
% 
if ifplot
[orisigpoints, bootsigpoints, clustersigpoints]=cluster_based_mcc_ttest(squeeze(betavalue)',0.05,1000);
% if nansum(orisigpoints)>0
%     plot(timelab,-orisigpoints*0.01,'color',[0.64,0.58,0.5]);
%     hold on
% end
% 
% if nansum(bootsigpoints)>0
%     plot(timelab,-bootsigpoints*0.005,'color',[0.96,0.64,0.38]);
%     hold on
% end
    if nansum(clustersigpoints)>0
        j=plot(timelab,clustersigpoints*0.005*(nplot-1.5),'color',setcolor,'linewidth',2);
set(get(get(j,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
ylabel(ylab,'FontSize',16);
xlabel('Time in ms from delivery of outcome','FontSize',16)
title(titlename);


