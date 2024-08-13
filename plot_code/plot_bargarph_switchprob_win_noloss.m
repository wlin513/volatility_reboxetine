switchprob_win_noloss=[squeeze(switchprob_win.(group{1})(v,:,:)),squeeze(switchprob_noloss.(group{1})(v,:,:))];
mean_switchprob=mean(switchprob_win_noloss);
sem_switchprob=std(switchprob_win_noloss)./sqrt(size(switchprob_win_noloss,1));
mean_switchprob=reshape(mean_switchprob,4,2);
sem_switchprob=reshape(sem_switchprob,4,2);
cond1={'win','noloss'};
f1=figure;
H=bar(mean_switchprob);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Switch probability','FontSize',16);
legend('win','noloss','AutoUpdate','off');
hold on
pos=1;
for ii=1:length(cond1)
for iii=1:length(blkname)
  pos=iii+0.14*(-1)^ii;
  plot([pos,pos],[mean_switchprob(iii,ii)-sem_switchprob(iii,ii),mean_switchprob(iii,ii)+sem_switchprob(iii,ii)],'-k','LineWidth',1)
end
end
hold off

% hold on
% for i=1:size(switchprob_win_noloss,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=switchprob_win_noloss(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% end
% end
% hold off
title([group{1},':visit',num2str(v)])
saveas(f1,[figdir,group{1},'_visit_',num2str(v),'_switchprob_win_noloss_result.png'])
%% plotting adaptation from stable to vol for other -vol stable or vol
switchprob_win_noloss=[switchprob_win.(group{1})(v,:,1)-switchprob_win.(group{1})(v,:,3);...
                                            switchprob_win.(group{1})(v,:,2)-switchprob_win.(group{1})(v,:,4);...
                                            switchprob_noloss.(group{1})(v,:,1)-switchprob_noloss.(group{1})(v,:,2);...
                                            switchprob_noloss.(group{1})(v,:,3)-switchprob_noloss.(group{1})(v,:,4)];
mean_switchprob=mean(switchprob_win_noloss,2);
sem_switchprob=std(switchprob_win_noloss,0,2)./sqrt(length(sublist.(group{1})));
mean_switchprob=reshape(mean_switchprob,2,2);
sem_switchprob=reshape(sem_switchprob,2,2);
cond1={'win','noloss'};
f1=figure;
H=bar(mean_switchprob);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
hold on
pos=1;
for ii=1:size(mean_switchprob,2)
for iii=1:size(mean_switchprob,1)
  pos=iii+0.15*(-1)^ii;
  plot([pos,pos],[mean_switchprob(iii,ii)-sem_switchprob(iii,ii),mean_switchprob(iii,ii)+sem_switchprob(iii,ii)],'-k','LineWidth',1)
end
end
hold off
set(gca,'XTickLabel',{'volatile','stable'},'FontSize',13);
xlabel('other schedule')
ylabel('Switch probability volatile vs stable','FontSize',16);
legend('win','noloss','AutoUpdate','off');
title([group{1},':visit',num2str(v)])
saveas(f1,[figdir,group{1},'_visit_',num2str(v),'_switchprob_vol_vs_sta_win_noloss_result.png'])