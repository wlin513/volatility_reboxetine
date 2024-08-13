bias_2lr_1b_bias_abs=squeeze(abs(bias_2lr_1b_bias.(group{1})(v,:,:)));
mean_bias_2lr_1b_bias_abs=mean(bias_2lr_1b_bias_abs,1);
sem_bias_2lr_1b_bias_abs=std(bias_2lr_1b_bias_abs,1)./sqrt(size(bias_2lr_1b_bias_abs,1));
up_bias_2lr_1b_bias=mean_bias_2lr_1b_bias_abs+sem_bias_2lr_1b_bias_abs;
low_bias_2lr_1b_bias=mean_bias_2lr_1b_bias_abs-sem_bias_2lr_1b_bias_abs;

mean_bias=[mean_bias_2lr_1b_bias_abs'];
up_bias_2lr1b_bias=[up_bias_2lr_1b_bias'];
low_bias_2lr1b_bias=[low_bias_2lr_1b_bias'];
f=figure;
H=bar(mean_bias);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('abs(bias)','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_bias_2lr1b_bias(iii),up_bias_2lr1b_bias(iii)],'-k','LineWidth',1)
end

% hold on
% for i=1:size(bias_2lr_1b_bias,1)
% for iii=1:length(blkname)
%   tmp=bias_2lr_1b_bias(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_bias_diff_result_2lr1b_bias.png'])