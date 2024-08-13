switchprob_neg_v1=[squeeze(switchprob_nowin.(group{1})(1,:,:)),squeeze(switchprob_loss.(group{1})(1,:,:))];
switchprob_neg_v2=[squeeze(switchprob_nowin.(group{1})(2,:,:)),squeeze(switchprob_loss.(group{1})(2,:,:))];
switchprob_neg=switchprob_neg_v2-switchprob_neg_v1;
mean_switchprob=mean(switchprob_neg);
sem_switchprob=std(switchprob_neg)./sqrt(size(switchprob_neg,1));
mean_switchprob=reshape(mean_switchprob,4,2);
sem_switchprob=reshape(sem_switchprob,4,2);
cond1={'nowin','loss'};
f1=figure;
H=bar(mean_switchprob);
H(1).FaceColor = 'green';
H(2).FaceColor = 'red';
set(gca,'XTickLabel',blkname,'FontSize',10);
ylabel('negative switch','FontSize',12);
legend('nowin','loss','AutoUpdate','off','Location','NorthWest');
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
    scatter(i*ones(1,length(sublist.(group{1})))-0.15,squeeze(switchprob_nowin.(group{1})(2,:,i)-switchprob_nowin.(group{1})(1,:,i)))
    hold on
    scatter(i*ones(1,length(sublist.(group{1})))+0.15,squeeze(switchprob_loss.(group{1})(2,:,i)-switchprob_loss.(group{1})(1,:,i)))
    hold on
    end
end
hold off

% hold on
% for i=1:size(switchprob_neg,1)
% for ii=1:length(cond1)
% for iii=1:length(blkname)
%   tmp=switchprob_neg(i,4*(ii-1)+iii);
%   pos=iii+0.14*(-1)^ii;
%   plot(pos,tmp,'xb')
% end
% endtitle([group{1},':visit',num2str(v)])
% end
% hold off
title([group{1},':visit2 vs visit1'])
saveas(f1,[figdir,group{1},'_visit2_vs visit1','_switchprob_neg_result.png'])