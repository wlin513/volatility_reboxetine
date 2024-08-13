win_vol_adpt_tr_diff.reboxetine=R_vol_adpt_othervol_model.win_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.win_vol_adpt_tr.reboxetine(1,:);
win_sta_adpt_tr_diff.reboxetine=R_vol_adpt_othervol_model.win_sta_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.win_sta_adpt_tr.reboxetine(1,:);
win_vol_adpt_tr_diff.placebo=R_vol_adpt_othervol_model.win_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.win_vol_adpt_tr.placebo(1,:);
win_sta_adpt_tr_diff.placebo=R_vol_adpt_othervol_model.win_sta_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.win_sta_adpt_tr.placebo(1,:);

means=[mean(win_vol_adpt_tr_diff.reboxetine),mean(win_vol_adpt_tr_diff.placebo);...
    mean(win_sta_adpt_tr_diff.reboxetine),mean(win_sta_adpt_tr_diff.placebo)];
sems=[std(win_vol_adpt_tr_diff.reboxetine)./sqrt(length(win_vol_adpt_tr_diff.reboxetine)),...
    std(win_vol_adpt_tr_diff.placebo)./sqrt(length(win_vol_adpt_tr_diff.placebo));...
    std(win_sta_adpt_tr_diff.reboxetine)./sqrt(length(win_vol_adpt_tr_diff.reboxetine)),...
    std(win_sta_adpt_tr_diff.placebo)./sqrt(length(win_vol_adpt_tr_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
%ylim([-0.1 0.1])
%yticks([-0.1:0.04:0.1])
set(gca,'XTickLabel',{'other(loss) volatile','other(loss) stable'},'FontSize',13);
ylabel('win learning rate volatility adaptation visit2 vs. visit1','FontSize',13);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');

hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pct_barplot_win_lr_adpt_othervol_diff_reb_pla.png'])
%%
loss_vol_adpt_tr_diff.reboxetine=R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(1,:);
loss_sta_adpt_tr_diff.reboxetine=R_vol_adpt_othervol_model.loss_sta_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_sta_adpt_tr.reboxetine(1,:);
loss_vol_adpt_tr_diff.placebo=R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(1,:);
loss_sta_adpt_tr_diff.placebo=R_vol_adpt_othervol_model.loss_sta_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_sta_adpt_tr.placebo(1,:);

means=[mean(loss_vol_adpt_tr_diff.reboxetine),mean(loss_vol_adpt_tr_diff.placebo);...
    mean(loss_sta_adpt_tr_diff.reboxetine),mean(loss_sta_adpt_tr_diff.placebo)];
sems=[std(loss_vol_adpt_tr_diff.reboxetine)./sqrt(length(loss_vol_adpt_tr_diff.reboxetine)),...
    std(loss_vol_adpt_tr_diff.placebo)./sqrt(length(loss_vol_adpt_tr_diff.placebo));...
    std(loss_sta_adpt_tr_diff.reboxetine)./sqrt(length(loss_vol_adpt_tr_diff.reboxetine)),...
    std(loss_sta_adpt_tr_diff.placebo)./sqrt(length(loss_vol_adpt_tr_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
%ylim([-0.1 0.1])
%yticks([-0.1:0.04:0.1])
set(gca,'XTickLabel',{'other(win) volatile','other(win) stable'},'FontSize',13);
ylabel('loss learning rate volatility adaptation visit2 vs. visit1','FontSize',13);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');

hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pct_barplot_loss_lr_adpt_othervol_diff_reb_pla.png'])