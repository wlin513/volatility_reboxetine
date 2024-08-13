alpha_1lr_1b_inv=squeeze(inv_logit(alpha_1lr_1b.(group{1})(v,:,:)));
mean_alpha_1lr_1b_inv=mean(alpha_1lr_1b_inv,1);
sem_alpha_1lr_1b_inv=std(alpha_1lr_1b_inv,1)./sqrt(size(alpha_1lr_1b_inv,1));
mean_alpha_1lr_1b=inv_logit(mean_alpha_1lr_1b_inv,1);
up_alpha_1lr_1b=inv_logit(mean_alpha_1lr_1b_inv+sem_alpha_1lr_1b_inv,1);
low_alpha_1lr_1b=inv_logit(mean_alpha_1lr_1b_inv-sem_alpha_1lr_1b_inv,1);

mean_alpha=[mean_alpha_1lr_1b'];
up_alpha_1lr1b=[up_alpha_1lr_1b'];
low_alpha_1lr1b=[low_alpha_1lr_1b'];
f=figure;
H=bar(mean_alpha);
H.FaceColor = 'blue';
set(gca,'XTickLabel',blkname,'FontSize',13);
ylabel('Learning Rate','FontSize',16);
hold on

pos=1;
for iii=1:length(blkname)
  pos=iii;
  plot([pos,pos],[low_alpha_1lr1b(iii),up_alpha_1lr1b(iii)],'-k','LineWidth',1)
end


% hold on
% for i=1:size(alpha_1lr_1b,1)
% for iii=1:length(blkname)
%   tmp=alpha_1lr_1b(i,iii);
%   pos=iii;
%   plot(pos,tmp,'xb')
% end
% end
hold off
title([group{1},':  visit',num2str(v)])
saveas(f,[figdir,group{1},'_visit',num2str(v),'_alpha_result_1lr1b.png'])