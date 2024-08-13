%% plotting adaptation from stable to vol for other -vol stable or vol
switchprob_win_noloss=[switchprob_win.(group{1})(2,:,1)-switchprob_win.(group{1})(2,:,3)-switchprob_win.(group{1})(1,:,1)+switchprob_win.(group{1})(1,:,3);...
                                            switchprob_win.(group{1})(2,:,2)-switchprob_win.(group{1})(2,:,4)-switchprob_win.(group{1})(1,:,2)+switchprob_win.(group{1})(1,:,4);...
                                            switchprob_noloss.(group{1})(2,:,1)-switchprob_noloss.(group{1})(2,:,2)-switchprob_noloss.(group{1})(1,:,1)+switchprob_noloss.(group{1})(1,:,2);...
                                            switchprob_noloss.(group{1})(2,:,3)-switchprob_noloss.(group{1})(2,:,4)-switchprob_noloss.(group{1})(1,:,3)+switchprob_noloss.(group{1})(1,:,4)];
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
title(['visit-2 vs visit1'])
title([group{1},': visit2 vs visit1'])
saveas(f1,[figdir,group{1},'diff_switchprob_vol_vs_sta_win_noloss_result.png'])