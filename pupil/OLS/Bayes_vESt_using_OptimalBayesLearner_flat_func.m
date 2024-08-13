clear all
clc
getfolders
parentdir='/home/wlin/Documents/2017_fMRI_pilot';


% subs={'4','5','6','7','8','9','10','11','12','14','15','16','17','18','19','20','21','22',...
%      '23','24','25','26','27','28','29','30','31','32'};

 subs={'25','26','27','28','29','30','31','32'};
 
 flattenpoints=[80,160,240,320];
 %blockorder : both_vol(1),win_vol(2),loss_vol(3),both_stable(4)
blockorders1=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,1,4,2;3,4,1,2];%blockorder for the not reversed version
blockorders2=[1,3,2,4; 1,2,3,4;4,3,2,1;4,2,3,1;3,1,4,2;3,4,1,2;2,1,4,3;2,4,1,3];%blockorder for the reversed version

 for sub=1:size(subs,2)
 %loads the behavior data
    name=subs{sub};
    beh_pth=[parentdir,'/behavior/data/matfiles/'];
    load([beh_pth,name,'.mat']);
    
    ver=data.ifreverse;
    blktype=data.blktype;
    %decide which version of counterblkorder to use
    if ver==1
        blkorder=blockorders1(blktype,:);
    else
        blkorder=blockorders2(blktype,:);
    end
    
    tmprewoutbayes=OptimalBayesLearner_flat_func(data.winpos);
    tmprewoutbayes.sub=name;
    tmplossoutbayes=OptimalBayesLearner_flat_func(data.losspos);
    tmprewoutbayes.sub=name;
    
    rewBayesout(sub)=tmprewoutbayes;
    lossBayesout(sub)=tmplossoutbayes;
   
%     f1=figure
%     plot(rewBayesout(sub).vEst,'g')
%     hold on
%     plot(lossBayesout(sub).vEst,'r')
%     title(['s',subs{sub},'',num2str(blkorder)])
%     saveas(f1,[figdir,'vEst/s',subs{sub},'_vEst',num2str(data.blktype),'_',num2str(data.ifreverse),'.png'])

    clear tmprewoutbayes tmplossoutbayes data
 end
    
 save([datadir,'vEST_k0_2000_noflatten.mat'])
 