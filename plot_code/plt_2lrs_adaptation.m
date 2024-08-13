function plt_2lrs_adaptation(win_lrs,loss_lrs,group,v,modelname,figdir,sublist)
%calculate the mean of alpha in logit space and inverse it back to the
%original space
win_lrs_inv=[inv_logit(win_lrs.(group{1})(v,:,1))-inv_logit(win_lrs.(group{1})(v,:,3));...
    inv_logit(win_lrs.(group{1})(v,:,2))-inv_logit(win_lrs.(group{1})(v,:,4))];
mean_win_lrs_inv=mean(win_lrs_inv,2);
sem_win_lrs_inv=std(win_lrs_inv,0,2)./sqrt(length(sublist.(group{1})));
up_win_lrs_inv=mean_win_lrs_inv+sem_win_lrs_inv;
low_win_lrs_inv=mean_win_lrs_inv-sem_win_lrs_inv;

loss_lrs_inv=[inv_logit(loss_lrs.(group{1})(v,:,1))-inv_logit(loss_lrs.(group{1})(v,:,2));...
    inv_logit(loss_lrs.(group{1})(v,:,3))-inv_logit(loss_lrs.(group{1})(v,:,4))];
mean_loss_lrs_inv=mean(loss_lrs_inv,2);
sem_loss_lrs_inv=std(loss_lrs_inv,0,2)./sqrt(length(sublist.(group{1})));
up_loss_lrs_inv=mean_loss_lrs_inv+sem_loss_lrs_inv;
low_loss_lrs_inv=mean_loss_lrs_inv-sem_loss_lrs_inv;

%% plot bargarph
mean_alpha=[mean_win_lrs_inv,mean_loss_lrs_inv];
up_alpha_2lr1b=[up_win_lrs_inv,up_loss_lrs_inv];
low_alpha_2lr1b=[low_win_lrs_inv,low_loss_lrs_inv];
cond={'alpha_rew','alpha_loss'};
f=figure;
H=bar(mean_alpha);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',{'volatile','stable'},'FontSize',13);
xlabel('other schedule')
ylabel('Learning Rate volatile vs. stable','FontSize',16);
legend('Win \alpha','Loss \alpha','AutoUpdate','off');
hold on

pos=1;
for ii=1:size(mean_alpha,2)
for iii=1:size(mean_alpha,1)
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
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_diff_vol_vs_sta_alpha_result_',modelname,'.png'])
