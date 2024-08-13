load('olsresult_otherstable.mat')
lossbeta_other_stable=lossbeta;
rewbeta_other_stable=rewbeta;
load('olsresult_othervol.mat')
lossbeta_other_vol=lossbeta;
rewbeta_other_vol=rewbeta;
loss_pupil_diff.reboxetine(:,1)=mean(lossbeta_other_vol.reboxetine(2,:,7,1001:3500),4)-mean(lossbeta_other_vol.reboxetine(1,:,7,1001:3500),4);
loss_pupil_diff.reboxetine(:,2)=mean(lossbeta_other_stable.reboxetine(2,:,7,1001:3500),4)-mean(lossbeta_other_stable.reboxetine(1,:,7,1001:3500),4);
loss_pupil_diff.placebo(:,1)=mean(lossbeta_other_vol.placebo(2,:,7,1001:3500),4)-mean(lossbeta_other_vol.placebo(1,:,7,1001:3500),4);
loss_pupil_diff.placebo(:,2)=mean(lossbeta_other_stable.placebo(2,:,7,1001:3500),4)-mean(lossbeta_other_stable.placebo(1,:,7,1001:3500),4);
rew_pupil_diff.reboxetine(:,1)=mean(rewbeta_other_vol.reboxetine(2,:,7,1001:3500),4)-mean(rewbeta_other_vol.reboxetine(1,:,7,1001:3500),4);
rew_pupil_diff.reboxetine(:,2)=mean(rewbeta_other_stable.reboxetine(2,:,7,1001:3500),4)-mean(rewbeta_other_stable.reboxetine(1,:,7,1001:3500),4);
rew_pupil_diff.placebo(:,1)=mean(rewbeta_other_vol.placebo(2,:,7,1001:3500),4)-mean(rewbeta_other_vol.placebo(1,:,7,1001:3500),4);
rew_pupil_diff.placebo(:,2)=mean(rewbeta_other_stable.placebo(2,:,7,1001:3500),4)-mean(rewbeta_other_stable.placebo(1,:,7,1001:3500),4);
means=[mean(loss_pupil_diff.reboxetine)',mean(loss_pupil_diff.placebo)'];
sems=[transpose(std(loss_pupil_diff.reboxetine)./length(loss_pupil_diff.reboxetine)),...
    transpose(std(loss_pupil_diff.placebo)./length(loss_pupil_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
set(gca,'XTickLabel',{'other(win) volatile','other(win) stable'},'FontSize',13);
ylabel({'selfPE * self volatility on loss pupil responses';' visit2 vs. visit1'},'FontSize',16);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');

hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pupil_barplot_loss_selfPEselfvol_othervol_diff_reb_pla.png'])

means=[mean(rew_pupil_diff.reboxetine)',mean(rew_pupil_diff.placebo)'];
sems=[transpose(std(rew_pupil_diff.reboxetine)./length(rew_pupil_diff.reboxetine)),...
    transpose(std(rew_pupil_diff.placebo)./length(rew_pupil_diff.placebo))];
up_pcts=means+sems;
low_pcts=means-sems;
f=figure('Position',[71,424,785,475]);
H=bar(means);
H(1).FaceColor = [0.188 0.663, 0.871];
H(2).FaceColor = [0.937,0.863, 0.020];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
set(gca,'XTickLabel',{'other(loss) volatile','other(loss) stable'},'FontSize',13);
ylabel({'selfPE * self volatility on win pupil responses';' visit2 vs. visit1'},'FontSize',16);
legend('reboxetine','placebo','Location','NorthEastOutside','AutoUpdate','off');

hold on

pos=1;
for ii=1:2
for iii=1:2
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_pcts(iii,ii),up_pcts(iii,ii)],'-k','LineWidth',1)
end
end
saveas(f,[figdir,'pupil_barplot_win_selfPEselfvol_othervol_diff_reb_pla.png'])

%%
groups=[ones(length(loss_pupil_diff.reboxetine),1);ones(length(loss_pupil_diff.placebo),1)*2];
table_pupil=array2table([[loss_pupil_diff.reboxetine;loss_pupil_diff.placebo],...
                       [rew_pupil_diff.reboxetine;rew_pupil_diff.placebo],groups],...
    'VariableNames',{'loss_pupil_diff_vol','loss_pupil_diff_sta','win_pupil_diff_vol','win_pupil_diff_sta','group'});
 writetable(table_pupil,[datadir,'pupil_betas'])
%% 
f=figure;
ppil=[loss_pupil_diff.reboxetine;loss_pupil_diff.placebo];
lr=[R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(1,:),...
    R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(1,:)];
s=scatter(ppil(:,1),lr,0.001,'.')
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
l=lsline;
l.LineWidth=2;
l.Color='red';
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
ppil=[loss_pupil_diff.reboxetine];
lr=[R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(1,:)];
%[r,p]=corr(lr(:,1),ppil')
l1=scatter(ppil(:,1),lr',30,[0.188 0.663, 0.871],'filled')

%
hold on
ppil=[loss_pupil_diff.placebo];
lr=[R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(1,:)];
%[r,p]=corr(lr(:,1),ppil')
l2=scatter(ppil(:,1),lr',30,[0.937,0.863, 0.020],'filled')
legend('reboxetine','placebo','AutoUpdate','off')
%
%[r,p]=corr(lr(:,1),ppil')

hold on
text(0.12, -0.5,'r=-0.30*','FontSize',13)
%text(1.8, -0.2,['r=',num2str(round(corr(lr(:,1),ppil'),2)),'*'],'FontSize',13)

xlabel({'selfPE * self volatility on loss pupil responses',' visit2 vs. visit1'})
ylabel({'loss learning rate adaptation',' visit2 vs. visit1'})
saveas(f,[figdir,'corr_loss_selfPEselfvol_othervol_diff_and_loss_lr_adpt.png'])