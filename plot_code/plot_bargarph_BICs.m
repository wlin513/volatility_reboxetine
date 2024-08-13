function plot_bargarph_BICs(modelnames,blkname,xlm,varargin)
nmodel=nargin-3;
for i=1:nmodel
    means(:,i)=mean(varargin{i}.session1.all.all);
    sems(:,i)=std(varargin{i}.session1.all.all)./sqrt(size(varargin{i}.session1.all.all,1));
end


f=figure;
%%
b=barh(means);
legend(modelnames,'Location','northeastoutside','AutoUpdate','off');

xlabel('BIC estimates','FontSize',16);
set(gca,'YTickLabel',blkname,'FontSize',10);
ytickangle(90);
xlim(xlm)
xpos=cell2mat(get(b,'XData'))'+[b.XOffset];
hold on
er=errorbar(means,xpos,sems,'horizontal');
for i=1:nmodel
er(i).LineStyle='none';
er(i).Color='k';
end
%%
saveas(f,[figdir,group{1},'visit',num2str(v),'_mean_BIC_in_allmodels.png'])