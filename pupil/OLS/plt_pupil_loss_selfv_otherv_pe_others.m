%% loss selfvol*others*PE
f=figure;
plot_ev(lossbeta.reboxetine(1,:,7,501:end),timelab,'','',[0.3 0.9, 1],true,2)
hold on
plot_ev(lossbeta.reboxetine(2,:,7,501:end),timelab,'Reboxetine','Coeffients for self-vol*loss-PE',[0 0.5, 1],true,1)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_selfv_otherv_PE_reboxetine_others.png'])

f=figure;
plot_ev(lossbeta.placebo(1,:,7,501:end),timelab,'','Parameter estimates',[1,0.95, 0.020],true,1)
hold on
plot_ev(lossbeta.placebo(2,:,7,501:end),timelab,'Placebo','Coeffients for self-vol*loss-PE',[0.937,0.702, 0.322],true,2)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_selfv_otherv_PE_placebo_others.png'])

f=figure;
plot_ev(lossbeta.reboxetine(2,:,7,501:end)-lossbeta.reboxetine(1,:,7,501:end),timelab,'','Parameter estimates visit2 vs. visit1',[0.188 0.663, 0.871],true,1)
hold on
plot_ev(lossbeta.placebo(2,:,7,501:end)-lossbeta.placebo(1,:,7,501:end),timelab,'visit2 vs. visit1: other-stable','Coeffients for self-vol*loss-PE',[0.937,0.863, 0.020],true,2)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('reboxetine','placebo')
saveas(f,[figdir,'loss_selfv_otherv_PE_reboxetine_v2minusv1_vs_placebo_others.png'])