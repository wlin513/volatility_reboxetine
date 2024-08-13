close all
clear;
clc;
cd('/home/wlin/Documents/2019_drug_study/code');
getfolders;
%%
addpath("plot_code/")
addpath("model-rescorla_wagner_2lr2b/")
%%
sub='2004';

datadir=[datadir,sub,'/'];

filename1=dir([datadir,'*_visit_1_blktype_*.txt']);
filename2=dir([datadir,'*_visit_2_blktype_*.txt']);

data(1)=read_txt([datadir,filename1.name]);
data(2)=read_txt([datadir,filename2.name]);

for i=1:2
tn=80;
blkn=4;
abandontn=10; % abandon first 10 trials when calculating posterior probability
abandontn_switch=0;
start=[0.5,0.5];

%all matrices for orders of 4 blocks: for 1 for both volatile, 2 for win
%volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];%blockorder for the not reversed version

blkname={'Both Volatile','Win Volatile','Loss Volatile','Both Stable'};
blktype=data(i).blktype;

%summaries the results for stay probolibities
for j=1:blkn
realblkname(j)=blkname(blockorders(blktype,j));
end

%% plot the schedules
for n=1:blkn
    est(((n-1)*tn+1):n*tn,:)=rescorla_wagner_2lr([data(i).winpos(((n-1)*tn+1):n*tn),data(i).losspos(((n-1)*tn+1):n*tn)],[0.3 0.3],start);
end
f3=figure;
%title(['sub:',num2str(data.subnum)]);
title(['sub:',sub,'   ','blk:',num2str(data(i).blktype)]);
hold on;
w1=[1 tn];w2=[(tn*2+1) tn*3];
h=[1 1];
a1=area(w1,h);a2=area(w2,h);
set(a1,'EdgeColor','none');set(a2,'EdgeColor','none');
set(a1,'FaceColor',[1 1 0.9]);set(a2,'FaceColor',[1 1 0.9]);
hold on;
    for u=1:blkn
        inferred_probs((((u-1)*tn+1):(u*tn)),:)=rescorla_wagner_2lr([data(i).winpos(((u-1)*tn+1):(u*tn)),data(i).losspos(((u-1)*tn+1):(u*tn))],[0.2 0.2],start);
    end
plot(inferred_probs(:,1),'g');
plot(inferred_probs(:,2),'r');
ylim([0 1])
ax = gca;
ax.XTick =[tn/2:tn:tn*4];
xticklabels(realblkname)
saveas(f3,[figdir,'schedule_sub_',sub,'.png']);
%saveas(f3,['schedule_sub',num2str(data.subnum),'.png']);
end