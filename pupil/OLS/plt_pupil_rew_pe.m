%%
f=figure;
plot_ev(rewbeta.reboxetine(1,:,2,501:end),timelab,'','',[0.3 0.9, 1],true,2)
hold on
plot_ev(rewbeta.reboxetine(2,:,2,501:end),timelab,'Reboxetine','Coeffients for win PE',[0 0.5, 1],true,1)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'rew_PE_reboxetine.png'])

f=figure;
plot_ev(rewbeta.placebo(1,:,2,501:end),timelab,'','',[1,0.95, 0.020],true,2)
hold on
plot_ev(rewbeta.placebo(2,:,2,501:end),timelab,'Placebo','Coeffients for win PE',[0.937,0.702, 0.322],true,1)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'rew_PE_placebo.png'])

f=figure;
plot_ev(rewbeta.reboxetine(2,:,2,501:end)-rewbeta.reboxetine(1,:,2,501:end),timelab,'','',[0.188 0.663, 0.871])
hold on
plot_ev(rewbeta.placebo(2,:,2,501:end)-rewbeta.placebo(1,:,2,501:end),timelab,'visit2 vs. visit1','Coeffients for win PE',[0.937,0.863, 0.020])
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('reboxetine','placebo')
saveas(f,[figdir,'rew_PE_reboxetine_v2minusv1_vs_placebo.png'])