beta_2lr_1b_bias_log=squeeze(log(beta_2lr_1b_bias.(group{1})(v,:,:)));
mean_beta_2lr_1b_bias_log=mean(beta_2lr_1b_bias_log,1);
sem_beta_2lr_1b_bias_log=std(beta_2lr_1b_bias_log,1)./sqrt(size(beta_2lr_1b_bias_log,1));
mean_beta_2lr_1b_bias=exp(mean_beta_2lr_1b_bias_log);
up_beta_2lr_1b_bias=exp(mean_beta_2lr_1b_bias_log+sem_beta_2lr_1b_bias_log);
low_beta_2lr_1b_bias=exp(mean_beta_2lr_1b_bias_log-sem_beta_2lr_1b_bias_log);

mean_beta=[mean_beta_2lr_1b_bias'];
up_beta_2lr1b_bias=[up_beta_2lr_1b_bias'];
low_beta_2lr1b_bias=[low_beta_2lr_1b_bias'];
f=figure;
H=bar(mean_beta);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Beta','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_beta_2lr1b_bias(iii),up_beta_2lr1b_bias(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(beta_2lr_1b_bias,1)
% for iii=1:length(blkname)
%   tmp=beta_2lr_1b_bias(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_beta_result_2lr1b_bias.png'])