switchprob_loss_psns_adpt_winv_diff.reboxetine=mean(switchprob_loss_psns.reboxetine(2,:,1)-switchprob_loss_psns.reboxetine(2,:,2),3)...
    -mean(switchprob_loss_psns.reboxetine(1,:,1)-switchprob_loss_psns.reboxetine(1,:,2),3);
switchprob_loss_psns_adpt_wins_diff.reboxetine=mean(switchprob_loss_psns.reboxetine(2,:,3)-switchprob_loss_psns.reboxetine(2,:,4),3)...
    -mean(switchprob_loss_psns.reboxetine(1,:,3)-switchprob_loss_psns.reboxetine(1,:,4),3);
switchprob_loss_psns_adpt_winv_diff.placebo=mean(switchprob_loss_psns.placebo(2,:,1)-switchprob_loss_psns.placebo(2,:,2),3)...
    -mean(switchprob_loss_psns.placebo(1,:,1)-switchprob_loss_psns.placebo(1,:,2),3);
switchprob_loss_psns_adpt_wins_diff.placebo=mean(switchprob_loss_psns.placebo(2,:,3)-switchprob_loss_psns.placebo(2,:,4),3)...
    -mean(switchprob_loss_psns.placebo(1,:,3)-switchprob_loss_psns.placebo(1,:,4),3);
means=[mean(switchprob_loss_psns_adpt_winv_diff.reboxetine),mean(switchprob_loss_psns_adpt_winv_diff.placebo);...
    mean(switchprob_loss_psns_adpt_wins_diff.reboxetine),mean(switchprob_loss_psns_adpt_wins_diff.placebo)];
sems=[std(switchprob_loss_psns_adpt_winv_diff.reboxetine)./sqrt(length(switchprob_loss_psns_adpt_winv_diff.reboxetine)),...
    std(switchprob_loss_psns_adpt_winv_diff.placebo)./sqrt(length(switchprob_loss_psns_adpt_winv_diff.placebo));...
    std(switchprob_loss_psns_adpt_wins_diff.reboxetine)./sqrt(length(switchprob_loss_psns_adpt_winv_diff.reboxetine)),...
    std(switchprob_loss_psns_adpt_wins_diff.placebo)./sqrt(length(switchprob_loss_psns_adpt_winv_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
ylim([-0.1 0.1])
yticks([-0.1:0.04:0.1])
set(gca,'XTickLabel',{'other(win) volatile','other(win) stable'},'FontSize',13);
ylabel('Ppsns loss volatility adaptation visit2 vs. visit1','FontSize',16);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');

hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pct_barplot_loss_psns_adpt_othervol_diff_reb_pla.png'])