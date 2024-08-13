switchprob_PNs_v1=[squeeze(1-switchprob_win.(group{1})(1,:,:)-switchprob_nowin.(group{1})(1,:,:)),squeeze(1-switchprob_noloss.(group{1})(1,:,:)-switchprob_loss.(group{1})(1,:,:))];
switchprob_PNs_v2=[squeeze(1-switchprob_win.(group{1})(2,:,:)-switchprob_nowin.(group{1})(2,:,:)),squeeze(1-switchprob_noloss.(group{1})(2,:,:)-switchprob_loss.(group{1})(2,:,:))];
switchprob_PNs=switchprob_PNs_v2-switchprob_PNs_v1;
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
legend('winPN','lossPN','AutoUpdate','off','Location','SouthWest');
hold on
pos=1;
for ii=1:length(cond1)
for iii=1:length(blkname)
  pos=iii+0.15*(-1)^ii;
  plot([pos,pos],[mean_switchprob(iii,ii)-sem_switchprob(iii,ii),mean_switchprob(iii,ii)+sem_switchprob(iii,ii)],'-k','LineWidth',1)
end
end

if ifplotscatter
    for i=1:4
    scatter(i*ones(1,length(sublist.(group{1})))-0.15,squeeze(1-switchprob_win.(group{1})(2,:,i)-switchprob_nowin.(group{1})(2,:,i)-(1-switchprob_win.(group{1})(1,:,i)-switchprob_nowin.(group{1})(1,:,i))))
    hold on
    scatter(i*ones(1,length(sublist.(group{1})))+0.15,squeeze(1-switchprob_noloss.(group{1})(2,:,i)-switchprob_loss.(group{1})(2,:,i)-(1-switchprob_noloss.(group{1})(1,:,i)-switchprob_loss.(group{1})(1,:,i))))
    hold on
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
title([group{1},':visit2 vs visit1'])
saveas(f1,[figdir,group{1},'_visit2_vs visit1','_switchprob_PNs_result.png'])