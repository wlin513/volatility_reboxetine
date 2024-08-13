function plt_2epsilon_dif(epsilon1,epsilon2,group,modelname,figdir)
blkname={'both volatile','win volatile','loss volatile','both stable'};
%calculate the mean of epsilon in logit space and storeerse it back to the
%original space
epsilon1_store=squeeze(epsilon1.(group{1})(2,:,:)-epsilon1.(group{1})(1,:,:));
mean_epsilon1=mean(epsilon1_store,1);
sem_epsilon1=std(epsilon1_store,1)./sqrt(size(epsilon1_store,1));
up_epsilon1=mean_epsilon1+sem_epsilon1;
low_epsilon1=mean_epsilon1-sem_epsilon1;

epsilon2_store=squeeze(epsilon2.(group{1})(2,:,:)-epsilon2.(group{1})(1,:,:));
mean_epsilon2=mean(epsilon2_store,1);
sem_epsilon2=std(epsilon2_store,1)./sqrt(size(epsilon2_store,1));
up_epsilon2=mean_epsilon2+sem_epsilon2;
low_epsilon2=mean_epsilon2-sem_epsilon2;
    

%% plot bargarph
mean_epsilon=[mean_epsilon1',mean_epsilon2'];
up_epsilon_2lr1b=[up_epsilon1',up_epsilon2'];
low_epsilon_2lr1b=[low_epsilon1',low_epsilon2'];
cond={'epsilon_rew','epsilon_loss'};
f=figure;
H=bar(mean_epsilon);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('epsilon','FontSize',16);
legend('Win \epsilon','Loss \epsilon','AutoUpdate','off');
hold on

pos=1;
for ii=1:length(cond)
for iii=1:length(blkname)
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[low_epsilon_2lr1b(iii,ii),up_epsilon_2lr1b(iii,ii)],'-k','LineWidth',1)
end
end


% epsilons_2lr_1b=[epsilon1,epsilon2];
% hold on
% for i=1:size(epsilons_2lr_1b,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=epsilons_2lr_1b(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% end
% end
hold off
title([group{1},':  epsilon visit2 - visit1'])
saveas(f,[figdir,group{1},'_epsilon_dif_result_',modelname,'.png'])
