pct_winchosen_mean_diff.reboxetine=mean(pct_winchosen.reboxetine(2,:,:),3)-mean(pct_winchosen.reboxetine(1,:,:),3);
pct_lossnotchosen_mean_diff.reboxetine=mean(pct_lossnotchosen.reboxetine(2,:,:),3)-mean(pct_lossnotchosen.reboxetine(1,:,:),3);
pct_winchosen_mean_diff.placebo=mean(pct_winchosen.placebo(2,:,:),3)-mean(pct_winchosen.placebo(1,:,:),3);
pct_lossnotchosen_mean_diff.placebo=mean(pct_lossnotchosen.placebo(2,:,:),3)-mean(pct_lossnotchosen.placebo(1,:,:),3);
means=[mean(pct_winchosen_mean_diff.reboxetine),mean(pct_winchosen_mean_diff.placebo);...
    mean(pct_lossnotchosen_mean_diff.reboxetine),mean(pct_lossnotchosen_mean_diff.placebo)];
sems=[std(pct_winchosen_mean_diff.reboxetine)./sqrt(length(pct_winchosen_mean_diff.reboxetine)),...
    std(pct_winchosen_mean_diff.placebo)./sqrt(length(pct_winchosen_mean_diff.placebo));...
    std(pct_lossnotchosen_mean_diff.reboxetine)./sqrt(length(pct_winchosen_mean_diff.reboxetine)),...
    std(pct_lossnotchosen_mean_diff.placebo)./sqrt(length(pct_winchosen_mean_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
ylim([-0.02,0.03])
yticks([-0.02:0.01:0.03])
set(gca,'XTickLabel',{'win','loss'},'FontSize',13);
ylabel('performance visit2 vs. visit1','FontSize',16);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');
hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pct_barplot_perf_diff_reb_pla.png'])