%calculate the mean of alpha in logit space and inverse it back to the
%original space
v=1;
alpha_rew_2lr_1b_bias_inv_v1=[inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,1))-inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,3));...
    inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,2))-inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,4))];
mean_alpha_rew_2lr1b_bias_inv_v1=mean(alpha_rew_2lr_1b_bias_inv_v1,2);
sem_alpha_rew_2lr1b_bias_inv_v1=std(alpha_rew_2lr_1b_bias_inv_v1,0,2)./sqrt(length(sublist.(group{1})));
up_alpha_rew_2lr1b_bias_inv_v1=mean_alpha_rew_2lr1b_bias_inv_v1+sem_alpha_rew_2lr1b_bias_inv_v1;
low_alpha_rew_2lr1b_bias_inv_v1=mean_alpha_rew_2lr1b_bias_inv_v1-sem_alpha_rew_2lr1b_bias_inv_v1;

v=2;
alpha_rew_2lr_1b_bias_inv_v2=[inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,1))-inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,3));...
    inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,2))-inv_logit(alpha_rew_2lr_1b_bias.(group{1})(v,:,4))];
mean_alpha_rew_2lr1b_bias_inv_v2=mean(alpha_rew_2lr_1b_bias_inv_v2,2);
sem_alpha_rew_2lr1b_bias_inv_v2=std(alpha_rew_2lr_1b_bias_inv_v2,0,2)./sqrt(length(sublist.(group{1})));
up_alpha_rew_2lr1b_bias_inv_v2=mean_alpha_rew_2lr1b_bias_inv_v2+sem_alpha_rew_2lr1b_bias_inv_v2;
low_alpha_rew_2lr1b_bias_inv_v2=mean_alpha_rew_2lr1b_bias_inv_v2-sem_alpha_rew_2lr1b_bias_inv_v2;

%% plot bargarph
mean_alpha=[mean_alpha_rew_2lr1b_bias_inv_v1,mean_alpha_rew_2lr1b_bias_inv_v2];
up_alpha_2lr1b_bias=[up_alpha_rew_2lr1b_bias_inv_v1,up_alpha_rew_2lr1b_bias_inv_v2];
low_alpha_2lr1b_bias=[low_alpha_rew_2lr1b_bias_inv_v1,low_alpha_rew_2lr1b_bias_inv_v2];
cond={'alpha_rew','alpha_rew'};
f=figure;
H=bar(mean_alpha);
H(1).FaceColor = 'green';
H(2).FaceColor = [0.561,0.737,0.580];
set(gca,'XTickLabel',{'volatile','stable'},'FontSize',13);
xlabel('other schedule')
ylabel('Learning Rate volatile vs. stable','FontSize',16);
legend('visit1','visit2','AutoUpdate','off');
hold on

pos=1;
for ii=1:size(mean_alpha,2)
for iii=1:size(mean_alpha,1)
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_alpha_2lr1b_bias(iii,ii),up_alpha_2lr1b_bias(iii,ii)],'-k','LineWidth',1)
end
end


% alphas_2lr_1b_bias=[alpha_rew_2lr_1b_bias,alpha_rew_2lr_1b_bias];
% hold on
% for i=1:size(alphas_2lr_1b_bias,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=alphas_2lr_1b_bias(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% end
% end
hold off
title([group{1},': win alpha'])
saveas(f,[figdir,group{1},'_cmp_visits_diff_vol_vs_sta_win_alpha_result_2lr1b_bias.png'])
