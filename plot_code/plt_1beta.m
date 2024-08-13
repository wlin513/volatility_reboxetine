function plt_1beta(beta,group,v,modelname,figdir)
blkname={'both volatile','win volatile','loss volatile','both stable'};
beta_log=squeeze(log(beta.(group{1})(v,:,:)));
mean_beta_log=mean(beta_log,1);
sem_beta_log=std(beta_log,1)./sqrt(size(beta_log,1));
mean_beta=exp(mean_beta_log);
up_beta=exp(mean_beta_log+sem_beta_log);
low_beta=exp(mean_beta_log-sem_beta_log);

mean_beta=[mean_beta'];
up_beta_2lr1b=[up_beta'];
low_beta_2lr1b=[low_beta'];
f=figure('Position',[329,219.5,642,420]);
H=bar(mean_beta);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Betas','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_beta_2lr1b(iii),up_beta_2lr1b(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(beta,1)
% for iii=1:length(blkname)
%   tmp=beta(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_beta_result_',modelname,'.png'])