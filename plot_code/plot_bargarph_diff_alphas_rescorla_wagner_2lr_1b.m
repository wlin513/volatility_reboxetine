%calculate the mean of alpha in logit space and inverse it back to the
%original space
alpha_rew_2lr_1b_inv=[inv_logit(alpha_rew_2lr_1b.(group{1})(v,:,1))-inv_logit(alpha_rew_2lr_1b.(group{1})(v,:,3));...
    inv_logit(alpha_rew_2lr_1b.(group{1})(v,:,2))-inv_logit(alpha_rew_2lr_1b.(group{1})(v,:,4))];
mean_alpha_rew_2lr1b_inv=mean(alpha_rew_2lr_1b_inv,2);
sem_alpha_rew_2lr1b_inv=std(alpha_rew_2lr_1b_inv,0,2)./sqrt(length(sublist.(group{1})));
up_alpha_rew_2lr1b_inv=mean_alpha_rew_2lr1b_inv+sem_alpha_rew_2lr1b_inv;
low_alpha_rew_2lr1b_inv=mean_alpha_rew_2lr1b_inv-sem_alpha_rew_2lr1b_inv;

alpha_loss_2lr_1b_inv=[inv_logit(alpha_loss_2lr_1b.(group{1})(v,:,1))-inv_logit(alpha_loss_2lr_1b.(group{1})(v,:,2));...
    inv_logit(alpha_loss_2lr_1b.(group{1})(v,:,3))-inv_logit(alpha_loss_2lr_1b.(group{1})(v,:,4))];
mean_alpha_loss_2lr1b_inv=mean(alpha_loss_2lr_1b_inv,2);
sem_alpha_loss_2lr1b_inv=std(alpha_loss_2lr_1b_inv,0,2)./sqrt(length(sublist.(group{1})));
up_alpha_loss_2lr1b_inv=mean_alpha_loss_2lr1b_inv+sem_alpha_loss_2lr1b_inv;
low_alpha_loss_2lr1b_inv=mean_alpha_loss_2lr1b_inv-sem_alpha_loss_2lr1b_inv;

%% plot bargarph
mean_alpha=[mean_alpha_rew_2lr1b_inv,mean_alpha_loss_2lr1b_inv];
up_alpha_2lr1b=[up_alpha_rew_2lr1b_inv,up_alpha_loss_2lr1b_inv];
low_alpha_2lr1b=[low_alpha_rew_2lr1b_inv,low_alpha_loss_2lr1b_inv];
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


% alphas_2lr_1b=[alpha_rew_2lr_1b,alpha_loss_2lr_1b];
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
saveas(f,[figdir,group{1},'_visit',num2str(v),'_diff_vol_vs_sta_alpha_result_2lr1b.png'])
