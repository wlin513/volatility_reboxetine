%%
f=figure;
plot_ev(lossbeta.reboxetine(1,:,2,501:end),timelab,'','',[0.3 0.9, 1],true,2)
hold on
plot_ev(lossbeta.reboxetine(2,:,2,501:end),timelab,'Reboxetine','Coeffients for loss PE',[0 0.5, 1],true,1)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_PE_reboxetine.png'])

f=figure;
plot_ev(lossbeta.placebo(1,:,2,501:end),timelab,'','',[1,0.95, 0.020],true,2)
hold on
plot_ev(lossbeta.placebo(2,:,2,501:end),timelab,'Placebo','Coeffients for loss PE',[0.937,0.702, 0.322],true,1)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_PE_placebo.png'])

f=figure;
plot_ev(lossbeta.reboxetine(2,:,2,501:end)-lossbeta.reboxetine(1,:,2,501:end),timelab,'','',[0.188 0.663, 0.871])
hold on
plot_ev(lossbeta.placebo(2,:,2,501:end)-lossbeta.placebo(1,:,2,501:end),timelab,'visit2 vs. visit1','Coeffients for loss PE',[0.937,0.863, 0.020])
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('reboxetine','placebo')
saveas(f,[figdir,'loss_PE_reboxetine_v2minusv1_vs_placebo.png'])