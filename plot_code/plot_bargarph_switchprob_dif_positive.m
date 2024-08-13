switchprob_pos_v1=[squeeze(1-switchprob_win.(group{1})(1,:,:)),squeeze(1-switchprob_noloss.(group{1})(1,:,:))];
switchprob_pos_v2=[squeeze(1-switchprob_win.(group{1})(2,:,:)),squeeze(1-switchprob_noloss.(group{1})(2,:,:))];
switchprob_pos=switchprob_pos_v2-switchprob_pos_v1;
mean_switchprob=mean(switchprob_pos);
sem_switchprob=std(switchprob_pos)./sqrt(size(switchprob_pos,1));
mean_switchprob=reshape(mean_switchprob,4,2);
sem_switchprob=reshape(sem_switchprob,4,2);
cond1={'win','noloss'};
f1=figure;
H=bar(mean_switchprob);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',blkname,'FontSize',10);
ylabel('positve stay','FontSize',12);
legend('win','noloss','AutoUpdate','off','Location','NorthWest');
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
    scatter(i*ones(1,length(sublist.(group{1})))-0.15,squeeze(1-switchprob_win.(group{1})(2,:,i)-1+switchprob_win.(group{1})(1,:,i)))
    hold on
    scatter(i*ones(1,length(sublist.(group{1})))+0.15,squeeze(1-switchprob_noloss.(group{1})(2,:,i)-1+switchprob_noloss.(group{1})(1,:,i)))
    hold on
    end
end
hold off

% hold on
% for i=1:size(switchprob_pos,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=switchprob_pos(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% endtitle([group{1},':visit',num2str(v)])
% end
% hold off
title([group{1},':visit2 vs visit1'])
saveas(f1,[figdir,group{1},'_visit2_vs visit1','_switchprob_pos_result.png'])