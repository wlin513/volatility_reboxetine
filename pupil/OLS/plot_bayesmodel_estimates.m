%%
clear
clc

getfolders

% load data
load([datadir,'vEST_k0_2000_noflatten.mat'])
%%
blktypes=[7,1,2,3,4,5,7,8];
vers=[1,2,2,2,2,2,2,2];
for i=1:length(rewBayesout)
    if vers(i)==1
        blkorder=blockorders1(blktypes(i),:);
    else
        blkorder=blockorders2(blktypes(i),:);
    end
    f1=figure
    plot(rewBayesout(i).vEst,'g')
    hold on
    plot(lossBayesout(i).vEst,'r')
    title(['block order: ',' ',num2str(blkorder)])
    saveas(f1,[figdir,'vEst/vEst_k0_20_flatten',num2str(blkorder),'.png'])
end
%%
%%
clear
clc

getfolders

% load data
load([datadir,'vEST_k0_2000.mat'])
blktypes=[7,1,2,3,4,5,7,8];
vers=[1,2,2,2,2,2,2,2];
for i=1:length(rewBayesout)
    if vers(i)==1
        blkorder=blockorders1(blktypes(i),:);
    else
        blkorder=blockorders2(blktypes(i),:);
    end
    f1=figure
    plot(rewBayesout(i).vEst,'g')
    hold on
    plot(lossBayesout(i).vEst,'r')
    title(['block order: ',' ',num2str(blkorder)])
    saveas(f1,[figdir,'vEst/vEst_k0_2000_noflatten_',num2str(blkorder),'.png'])
end