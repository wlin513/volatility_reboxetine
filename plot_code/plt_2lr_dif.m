function plt_2lr_dif(win_lrs,loss_lrs,group,modelname,figdir)
blkname={'both volatile','win volatile','loss volatile','both stable'};
%calculate the mean of alpha in logit space and inverse it back to the
%original space
win_lrs_inv=squeeze(inv_logit(win_lrs.(group{1})(2,:,:))-inv_logit(win_lrs.(group{1})(1,:,:)));
mean_win_lrs_inv=mean(win_lrs_inv,1);
sem_win_lrs_inv=std(win_lrs_inv,1)./sqrt(size(win_lrs_inv,1));
mean_win_lrs=mean_win_lrs_inv;
up_win_lrs=mean_win_lrs_inv+sem_win_lrs_inv;
low_win_lrs=mean_win_lrs_inv-sem_win_lrs_inv;

loss_lrs_inv=squeeze(inv_logit(loss_lrs.(group{1})(2,:,:))-inv_logit(loss_lrs.(group{1})(1,:,:)));
mean_loss_lrs_inv=mean(loss_lrs_inv,1);
sem_loss_lrs_inv=std(loss_lrs_inv,1)./sqrt(size(loss_lrs_inv,1));
mean_loss_lrs=mean_loss_lrs_inv;
up_loss_lrs=mean_loss_lrs_inv+sem_loss_lrs_inv;
low_loss_lrs=mean_loss_lrs_inv-sem_loss_lrs_inv;
    

%% plot bargarph
mean_alpha=[mean_win_lrs',mean_loss_lrs'];
up_alpha_2lr1b=[up_win_lrs',up_loss_lrs'];
low_alpha_2lr1b=[low_win_lrs',low_loss_lrs'];
cond={'alpha_rew','alpha_loss'};
f=figure;
H=bar(mean_alpha);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
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
% for i=1:size(alphas_2lr_1b,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=alphas_2lr_1b(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% end
% end
hold off
title([group{1},':  alpha visit2 - visit1'])
saveas(f,[figdir,group{1},'_alpha_dif_result_',modelname,'.png'])
