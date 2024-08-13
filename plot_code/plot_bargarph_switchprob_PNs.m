switchprob_PNs=[squeeze(1-switchprob_win.(group{1})(v,:,:)-switchprob_nowin.(group{1})(v,:,:)),squeeze(1-switchprob_noloss.(group{1})(v,:,:)-switchprob_loss.(group{1})(v,:,:))];
mean_switchprob=mean(switchprob_PNs);
sem_switchprob=std(switchprob_PNs)./sqrt(size(switchprob_PNs,1));
mean_switchprob=reshape(mean_switchprob,4,2);
sem_switchprob=reshape(sem_switchprob,4,2);
cond1={'win','loss'};
f1=figure;
H=bar(mean_switchprob);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',blkname,'FontSize',10);
ylabel('positve stay vs. negative switch probability','FontSize',12);
legend('winPN','lossPN','AutoUpdate','off','Location','NorthWest');
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
% for i=1:size(switchprob_PNs,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=switchprob_PNs(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% endtitle([group{1},':visit',num2str(v)])
% end
% hold off
title([group{1},':visit',num2str(v)])
saveas(f1,[figdir,group{1},'_visit_',num2str(v),'_switchprob_PNs_result.png'])
%% plotting adaptation from stable to vol for other -vol stable or vol
switchprob_PNs=[switchprob_winPN.(group{1})(v,:,1)-switchprob_winPN.(group{1})(v,:,3);...
                                            switchprob_winPN.(group{1})(v,:,2)-switchprob_winPN.(group{1})(v,:,4);...
                                            switchprob_lossPN.(group{1})(v,:,1)-switchprob_lossPN.(group{1})(v,:,2);...
                                            switchprob_lossPN.(group{1})(v,:,3)-switchprob_lossPN.(group{1})(v,:,4)];
mean_switchprob=mean(switchprob_PNs,2);
sem_switchprob=std(switchprob_PNs,0,2)./sqrt(length(sublist.(group{1})));
mean_switchprob=reshape(mean_switchprob,2,2);
sem_switchprob=reshape(sem_switchprob,2,2);
cond1={'win','loss'};
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
legend('winPN','lossPN','AutoUpdate','off','Location','SouthWest');
title([group{1},':visit',num2str(v)])
saveas(f1,[figdir,group{1},'_visit_',num2str(v),'_switchprob_vol_vs_sta_PNs_result.png'])