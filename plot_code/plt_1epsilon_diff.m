function plt_1epsilon_diff(epsilon,group,modelname,figdir)
blkname={'both volatile','win volatile','loss volatile','both stable'};
epsilon_store=squeeze(epsilon.(group{1})(2,:,:)-epsilon.(group{1})(1,:,:));
mean_epsilon=mean(epsilon_store,1);
sem_epsilon=std(epsilon_store,1)./sqrt(size(epsilon_store,1));
up_epsilon=mean_epsilon+sem_epsilon;
low_epsilon=mean_epsilon-sem_epsilon;

mean_epsilon=[mean_epsilon'];
up_epsilon_2lr1b=[up_epsilon'];
low_epsilon_2lr1b=[low_epsilon'];
f=figure;
H=bar(mean_epsilon);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('epsilon','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_epsilon_2lr1b(iii),up_epsilon_2lr1b(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(epsilon,1)
% for iii=1:length(blkname)
%   tmp=epsilon(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  epsilon visit2 - visit1'])
saveas(f,[figdir,group{1},'_epsilon_diff_result_',modelname,'.png'])