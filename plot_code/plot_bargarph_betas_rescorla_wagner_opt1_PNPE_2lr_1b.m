beta_PNPE_2lr_1b_opt1_log=squeeze(log(beta_PNPE_2lr_1b_opt1.(group{1})(v,:,:)));
mean_beta_PNPE_2lr_1b_opt1_log=mean(beta_PNPE_2lr_1b_opt1_log,1);
sem_beta_PNPE_2lr_1b_opt1_log=std(beta_PNPE_2lr_1b_opt1_log,1)./sqrt(size(beta_PNPE_2lr_1b_opt1_log,1));
mean_beta_PNPE_2lr_1b_opt1=exp(mean_beta_PNPE_2lr_1b_opt1_log);
up_beta_PNPE_2lr_1b_opt1=exp(mean_beta_PNPE_2lr_1b_opt1_log+sem_beta_PNPE_2lr_1b_opt1_log);
low_beta_PNPE_2lr_1b_opt1=exp(mean_beta_PNPE_2lr_1b_opt1_log-sem_beta_PNPE_2lr_1b_opt1_log);

mean_beta=[mean_beta_PNPE_2lr_1b_opt1'];
up_beta_2lr1b=[up_beta_PNPE_2lr_1b_opt1'];
low_beta_2lr1b=[low_beta_PNPE_2lr_1b_opt1'];
f=figure;
H=bar(mean_beta);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Beta','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_beta_2lr1b(iii),up_beta_2lr1b(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(beta_PNPE_2lr_1b_opt1,1)
% for iii=1:length(blkname)
%   tmp=beta_PNPE_2lr_1b_opt1(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_beta_result_PNPE_2lr1b_opt1.png'])