
mean_stk=squeeze(mean(stk_2lr_1b_choice_stickness.(group{1})(2,:,:)-stk_2lr_1b_choice_stickness.(group{1})(1,:,:)));
up_stk_2lr1b_choice_stickness=mean_stk'+std(squeeze(stk_2lr_1b_choice_stickness.(group{1})(2,:,:)-stk_2lr_1b_choice_stickness.(group{1})(1,:,:)))./sqrt(length(stk_2lr_1b_choice_stickness.(group{1})(v,:,:)));
low_stk_2lr1b_choice_stickness=mean_stk'-std(squeeze(stk_2lr_1b_choice_stickness.(group{1})(2,:,:)-stk_2lr_1b_choice_stickness.(group{1})(1,:,:)))./sqrt(length(stk_2lr_1b_choice_stickness.(group{1})(v,:,:)));
f=figure;
H=bar(mean_stk);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('stk','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_stk_2lr1b_choice_stickness(iii),up_stk_2lr1b_choice_stickness(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(stk_2lr_1b_choice_stickness,1)
% for iii=1:length(blkname)
%   tmp=stk_2lr_1b_choice_stickness(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},': stk visit2 - visit1'])
saveas(f,[figdir,group{1},'_stk_diff_result_2lr1b_choice_stickness.png'])