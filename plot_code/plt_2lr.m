function plt_2lr(win_lrs,loss_lrs,group,v,modelname,figdir)
blkname={'both volatile','win volatile','loss volatile','both stable'};
%calculate the mean of alpha in logit space and inverse it back to the
%original space
win_lrs_inv=squeeze(inv_logit(win_lrs.(group{1})(v,:,:)));
mean_win_lrs_inv=nanmean(win_lrs_inv,1);
sem_win_lrs_inv=nanstd(win_lrs_inv,1)./sqrt(sum(~isnan(win_lrs_inv),1));
mean_win_lrs=inv_logit(mean_win_lrs_inv,1);
up_win_lrs=inv_logit(mean_win_lrs_inv+sem_win_lrs_inv,1);
low_win_lrs=inv_logit(mean_win_lrs_inv-sem_win_lrs_inv,1);

loss_lrs_inv=squeeze(inv_logit(loss_lrs.(group{1})(v,:,:)));
mean_loss_lrs_inv=nanmean(loss_lrs_inv,1);
sem_loss_lrs_inv=nanstd(loss_lrs_inv,1)./sqrt(sum(~isnan(loss_lrs_inv),1));
mean_loss_lrs=inv_logit(mean_loss_lrs_inv,1);
up_loss_lrs=inv_logit(mean_loss_lrs_inv+sem_loss_lrs_inv,1);
low_loss_lrs=inv_logit(mean_loss_lrs_inv-sem_loss_lrs_inv,1);
    

%% plot bargarph
mean_alpha=[mean_win_lrs',mean_loss_lrs'];
up_alpha_2lr1b=[up_win_lrs',up_loss_lrs'];
low_alpha_2lr1b=[low_win_lrs',low_loss_lrs'];
cond={'alpha_win','alpha_loss'};
f=figure('Position',[71,424,785,475]);
H=bar(mean_alpha);
H(1).FaceColor = [0.702 0.847, 0.38];
H(2).FaceColor = [0.945,0.419,0.435];
H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Learning Rate','FontSize',16);
legend('Win \alpha','Loss \alpha','AutoUpdate','off');
hold on

pos=1;
for ii=1:length(cond)
for iii=1:length(blkname)
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_alpha_2lr1b(iii,ii),up_alpha_2lr1b(iii,ii)],'-k','LineWidth',1)
end
end


% alphas_2lr_1b=[win_lrs,loss_lrs];
% hold on
% for i=1:sum(~isnan(alphas_2lr_1b,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=alphas_2lr_1b(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_alphas_model_',modelname,'.png'])
