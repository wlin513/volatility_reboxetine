function plot_ev_block(betavalue,timelab,titlename,setcolor,markplace)
% output: a plot with ev beta means, std error and significant points marked
% input: betavalue,timelab,titlename,setcolor,markplace(2 for ylim max,1 for ylim min)
%
mean1=squeeze(mean(betavalue))';
sem=squeeze(std(betavalue)./sqrt(size(betavalue,1)));
sem=sem';
jbfill(timelab,mean1+sem,mean1-sem,setcolor,setcolor,0,0.1);
hold on
plot(timelab,mean1,'color',setcolor);
hold on
yl=ylim;
[orisigpoints, bootsigpoints, clustersigpoints]=cluster_based_mcc_ttest(squeeze(betavalue)',0.05,1000);

if nansum(clustersigpoints)>0
    plot(timelab,clustersigpoints*(yl(markplace)-(-1)^markplace*0.01),setcolor,'linewidth',2);
    hold on
end

xlim([min(timelab) max(timelab)]);
ylabel('parameter estimates');
title(titlename);


