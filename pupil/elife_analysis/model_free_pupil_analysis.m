%% analyse pupil timecourse model free based on outcome delivered
cd('D:/2019_drug_study/code/pupil/elife_analysis');
clear all
getfolders
clc
ascdatadir=[datadir,'ascfiles/'];
abandontn=10;%abandon first 10 trails for each block 
nblk=4;
tn=80;
%matrices for orders of 4 blocks: for 1 for both volatile, 2 for win volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];
%% load subjects and group information
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
sublist.placebo=sublist.all(Ques.treatment==2);
excludesub=[1019,2013,1011];%2013 a lot of nans in the results;1011 for nans for loss in win vol block
for group={'reboxetine','placebo'}
    for i=1:length(excludesub)
        if isempty(find(sublist.(group{1})==excludesub(i)))
            continue
        else
              excnum(i)=find(sublist.(group{1})==excludesub(i));
        end
    sublist.(group{1})(excnum(i))=[];
    end
end
%%
for group={'reboxetine','placebo'}
    ssublist=sublist.(group{1});
    for sub=1:size(ssublist,1)    
        name=num2str(ssublist(sub));
        sdatadir=[datadir,name,'/'];
        for visit=1:2
                %% load pupil data 
                file_name=['out_',name,'_visit_',num2str(visit),'_not_normalised.mat'];
                load([ascdatadir,file_name])
                load([ascdatadir,name,'_visit_', num2str(visit),'_pupil_exclude_trials.mat']);
                rew_include=1-rew_exclude'; % to exclude trials where more than 50% of the data were missing
                pun_include=1-pun_exclude';
                rew_data=out.rew_outcome;
                pun_data=out.pun_outcome;
                %loads the behavior data
                filename=dir([sdatadir,'*_visit_',num2str(visit),'_blktype_*.txt']);
                data=read_txt([sdatadir,filename.name]);
                ifchosen_pun=data.losschosen;
                ifchosen_rew=data.winchosen;
                rew_order=data.order;%1 for win first 0 for loss first
                pun_order=1-data.order;%1 for loss first 0 for win first
                %% baseline corrections 
                rewbaselines=mean(rew_data(:,501:1000),2); %every 500 represents 1 seconds
                rew_data=rew_data-rewbaselines;
                punbaselines=mean(pun_data(:,501:1000),2); 
                pun_data=pun_data-punbaselines;
                
                %% exclude trails and then calculate mean across trials for each trial type
                includetrials=repmat([false(abandontn,1);true(tn-abandontn,1)],nblk,1);
                rewnantrials=~isnan(mean(rew_data,2));
                punnantrials=~isnan(mean(pun_data,2));
                summary_nantrials.(group{1})(visit,sub,:)=[sum(isnan(mean(rew_data,2))),sum(isnan(mean(pun_data,2)))];
                % generate block index so that every subject's data is organised in the order of 1 both vol,2 win vol, 3 loss vol, 4 both stable
                blkorder=blockorders(data.blktype,:);
                for blk=1:nblk   
                        blkindex=false(size(rew_data,1),1);
                        blkloc=find(blkorder==blk);
                        blkindex((blkloc-1)*tn+1:blkloc*tn)=true;
                        rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                        rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                        means_rew_delivered.(group{1})(visit,sub,blk,:)=mean(rew_delivered,1); %take means for each timepoint
                        means_rew_undelivered.(group{1})(visit,sub,blk,:)=mean(rew_undelivered,1);
                        pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                        pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                        means_pun_delivered.(group{1})(visit,sub,blk,:)=mean(pun_delivered,1);
                        means_pun_undelivered.(group{1})(visit,sub,blk,:)=mean(pun_undelivered,1);                       

                        rew_bslines=rewbaselines(blkindex&includetrials&rewnantrials&rew_include);
                        means_rew_bslines.(group{1})(visit,sub,blk)=mean(rew_bslines,1);
                        pun_bslines=punbaselines(blkindex&includetrials&punnantrials&pun_include);
                        means_pun_bslines.(group{1})(visit,sub,blk)=mean(pun_bslines,1);
                end
                %%mean for each conditions
                %mean for win vol blocks
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==1|blkorder==2);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                delta_win_vol.(group{1})(visit,sub,:)=mean(rew_delivered,1)-mean(rew_undelivered,1);
               %mean for win stable blocks
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==3|blkorder==4);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                delta_win_stable.(group{1})(visit,sub,:)=mean(rew_delivered,1)-mean(rew_undelivered,1);
                
               %mean for loss vol blocks
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==1|blkorder==3);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                delta_loss_vol.(group{1})(visit,sub,:)=mean(pun_delivered,1)-mean(pun_undelivered,1);
                %mean for loss stable blocks
                blkindex=false(size(pun_data,1),1);
                blkloc=find(blkorder==2|blkorder==4);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                delta_loss_stable.(group{1})(visit,sub,:)=mean(pun_delivered,1)-mean(pun_undelivered,1);  
                
                %mean across winvol and lossvol blocks(the semi volatile blocks)
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==2|blkorder==3);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                delta_win_semi.(group{1})(visit,sub,:)=mean(rew_delivered,1)-mean(rew_undelivered,1);
                means_rew_delivered_semi.(group{1})(visit,sub,:)=mean(rew_delivered,1);
                means_rew_undelivered_semi.(group{1})(visit,sub,:)=mean(rew_undelivered,1);
                
                pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                delta_loss_semi.(group{1})(visit,sub,:)=mean(pun_delivered,1)-mean(pun_undelivered,1);
                means_pun_delivered_semi.(group{1})(visit,sub,:)=mean(pun_delivered,1);
                means_pun_undelivered_semi.(group{1})(visit,sub,:)=mean(pun_undelivered,1);
               %mean across both vol and both stable blocks(the balance volatile blocks)
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==1|blkorder==4);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                delta_win_bal.(group{1})(visit,sub,:)=mean(rew_delivered,1)-mean(rew_undelivered,1);
                means_rew_delivered_bal.(group{1})(visit,sub,:)=mean(rew_delivered,1);
                means_rew_undelivered_bal.(group{1})(visit,sub,:)=mean(rew_undelivered,1);
                
                pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                delta_loss_bal.(group{1})(visit,sub,:)=mean(pun_delivered,1)-mean(pun_undelivered,1);
                means_pun_delivered_bal.(group{1})(visit,sub,:)=mean(pun_delivered,1);
                means_pun_undelivered_bal.(group{1})(visit,sub,:)=mean(pun_undelivered,1);
        end
    end
    delta_win.(group{1})=means_rew_delivered.(group{1})-means_rew_undelivered.(group{1}); %pupil response to rewards when chosen vs unchosen rewarded
    delta_loss.(group{1})=means_pun_delivered.(group{1})-means_pun_undelivered.(group{1});
    diff_delta_win.(group{1})=delta_win.(group{1})(2,:,:,:)-delta_win.(group{1})(1,:,:,:);
    diff_delta_loss.(group{1})=delta_loss.(group{1})(2,:,:,:)-delta_loss.(group{1})(1,:,:,:);
    % visit 2 - visit 1
    diff_delta_win_vol.(group{1})=delta_win_vol.(group{1})(2,:,:)-delta_win_vol.(group{1})(1,:,:);
    diff_delta_win_stable.(group{1})=delta_win_stable.(group{1})(2,:,:)-delta_win_stable.(group{1})(1,:,:);
    diff_delta_loss_vol.(group{1})=delta_loss_vol.(group{1})(2,:,:)-delta_loss_vol.(group{1})(1,:,:);
    diff_delta_loss_stable.(group{1})=delta_loss_stable.(group{1})(2,:,:)-delta_loss_stable.(group{1})(1,:,:);  
    
    diff_rew_delivered.(group{1})=means_rew_delivered.(group{1})(2,:,:,:)-means_rew_delivered.(group{1})(1,:,:,:);  
    diff_rew_undelivered.(group{1})=means_rew_undelivered.(group{1})(2,:,:,:)-means_rew_undelivered.(group{1})(1,:,:,:);
    diff_pun_delivered.(group{1})=means_pun_delivered.(group{1})(2,:,:,:)-means_pun_delivered.(group{1})(1,:,:,:);  
    diff_pun_undelivered.(group{1})=means_pun_undelivered.(group{1})(2,:,:,:)-means_pun_undelivered.(group{1})(1,:,:,:);
   
    diff_delta_win_semi.(group{1})=delta_win_semi.(group{1})(2,:,:)-delta_win_semi.(group{1})(1,:,:);
    diff_means_rew_delivered_semi.(group{1})=means_rew_delivered_semi.(group{1})(2,:,:)-means_rew_delivered_semi.(group{1})(1,:,:);
    diff_means_rew_undelivered_semi.(group{1})=means_rew_undelivered_semi.(group{1})(2,:,:)-means_rew_undelivered_semi.(group{1})(1,:,:);
    diff_delta_loss_semi.(group{1})=delta_loss_semi.(group{1})(2,:,:)-delta_loss_semi.(group{1})(1,:,:);
    diff_means_pun_delivered_semi.(group{1})=means_pun_delivered_semi.(group{1})(2,:,:)-means_pun_delivered_semi.(group{1})(1,:,:);
    diff_means_pun_undelivered_semi.(group{1})=means_pun_undelivered_semi.(group{1})(2,:,:)-means_pun_undelivered_semi.(group{1})(1,:,:);
       
    diff_delta_win_bal.(group{1})=delta_win_bal.(group{1})(2,:,:)-delta_win_bal.(group{1})(1,:,:);
    diff_means_rew_delivered_bal.(group{1})=means_rew_delivered_bal.(group{1})(2,:,:)-means_rew_delivered_bal.(group{1})(1,:,:);
    diff_means_rew_undelivered_bal.(group{1})=means_rew_undelivered_bal.(group{1})(2,:,:)-means_rew_undelivered_bal.(group{1})(1,:,:);
    diff_delta_loss_bal.(group{1})=delta_loss_bal.(group{1})(2,:,:)-delta_loss_bal.(group{1})(1,:,:);
    diff_means_pun_delivered_bal.(group{1})=means_pun_delivered_bal.(group{1})(2,:,:)-means_pun_delivered_bal.(group{1})(1,:,:);
    diff_means_pun_undelivered_bal.(group{1})=means_pun_undelivered_bal.(group{1})(2,:,:)-means_pun_undelivered_bal.(group{1})(1,:,:);
    %% now take the means and standard errors of each timecourse
       %
        group_means_rew_delivered.(group{1})=squeeze(mean(means_rew_delivered.(group{1}),2));
        group_means_rew_undelivered.(group{1})=squeeze(mean(means_rew_undelivered.(group{1}),2));
        group_means_pun_delivered.(group{1})=squeeze(mean(means_pun_delivered.(group{1}),2));
        group_means_pun_undelivered.(group{1})=squeeze(mean(means_pun_undelivered.(group{1}),2));

        group_stes_rew_delivered.(group{1})=squeeze(std(means_rew_delivered.(group{1}),0,2))./sqrt(size(means_rew_delivered.(group{1}),2));
        group_stes_rew_undelivered.(group{1})=squeeze(std(means_rew_undelivered.(group{1}),0,2))./sqrt(size(means_rew_delivered.(group{1}),2));
        group_stes_pun_delivered.(group{1})=squeeze(std(means_pun_delivered.(group{1}),0,2))./sqrt(size(means_rew_delivered.(group{1}),2));
        group_stes_pun_undelivered.(group{1})=squeeze(std(means_pun_undelivered.(group{1}),0,2))./sqrt(size(means_rew_delivered.(group{1}),2));
         
        group_means_diff_rew_delivered.(group{1})=squeeze(mean(diff_rew_delivered.(group{1}),2));
        group_means_diff_rew_undelivered.(group{1})=squeeze(mean(diff_rew_undelivered.(group{1}),2));
        group_means_diff_pun_delivered.(group{1})=squeeze(mean(diff_pun_delivered.(group{1}),2));
        group_means_diff_pun_undelivered.(group{1})=squeeze(mean(diff_pun_undelivered.(group{1}),2));
        
        group_stes_diff_rew_delivered.(group{1})=squeeze(std(diff_rew_delivered.(group{1}),0,2))./sqrt(size(diff_rew_delivered.(group{1}),2));
        group_stes_diff_rew_undelivered.(group{1})=squeeze(std(diff_rew_undelivered.(group{1}),0,2))./sqrt(size(diff_rew_undelivered.(group{1}),2));
        group_stes_diff_pun_delivered.(group{1})=squeeze(std(diff_pun_delivered.(group{1}),0,2))./sqrt(size(diff_pun_delivered.(group{1}),2));
        group_stes_diff_pun_undelivered.(group{1})=squeeze(std(diff_pun_undelivered.(group{1}),0,2))./sqrt(size(diff_pun_undelivered.(group{1}),2));
  
        %
        group_means_delta_win_vol.(group{1})=squeeze(mean(delta_win_vol.(group{1}),2));
        group_stes_delta_win_vol.(group{1})=squeeze(std(delta_win_vol.(group{1}),0,2))./sqrt(size(delta_win_vol.(group{1}),2));
        group_means_delta_win_stable.(group{1})=squeeze(mean(delta_win_stable.(group{1}),2));
        group_stes_delta_win_stable.(group{1})=squeeze(std(delta_win_stable.(group{1}),0,2))./sqrt(size(delta_win_stable.(group{1}),2));
        group_means_delta_loss_vol.(group{1})=squeeze(mean(delta_loss_vol.(group{1}),2));
        group_stes_delta_loss_vol.(group{1})=squeeze(std(delta_loss_vol.(group{1}),0,2))./sqrt(size(delta_loss_vol.(group{1}),2));
        group_means_delta_loss_stable.(group{1})=squeeze(mean(delta_loss_stable.(group{1}),2));
        group_stes_delta_loss_stable.(group{1})=squeeze(std(delta_loss_stable.(group{1}),0,2))./sqrt(size(delta_loss_stable.(group{1}),2));
        
        group_means_diff_delta_win_vol.(group{1})=transpose(squeeze(mean(diff_delta_win_vol.(group{1}),2)));
        group_stes_diff_delta_win_vol.(group{1})=transpose(squeeze(std(diff_delta_win_vol.(group{1}),0,2))./sqrt(size(diff_delta_win_vol.(group{1}),2)));
        group_means_diff_delta_win_stable.(group{1})=transpose(squeeze(mean(diff_delta_win_stable.(group{1}),2)));
        group_stes_diff_delta_win_stable.(group{1})=transpose(squeeze(std(diff_delta_win_stable.(group{1}),0,2))./sqrt(size(diff_delta_win_stable.(group{1}),2)));
        group_means_diff_delta_loss_vol.(group{1})=transpose(squeeze(mean(diff_delta_loss_vol.(group{1}),2)));
        group_stes_diff_delta_loss_vol.(group{1})=transpose(squeeze(std(diff_delta_loss_vol.(group{1}),0,2))./sqrt(size(diff_delta_loss_vol.(group{1}),2)));
        group_means_diff_delta_loss_stable.(group{1})=transpose(squeeze(mean(diff_delta_loss_stable.(group{1}),2)));
        group_stes_diff_delta_loss_stable.(group{1})=transpose(squeeze(std(diff_delta_loss_stable.(group{1}),0,2))./sqrt(size(diff_delta_loss_stable.(group{1}),2)));

         %
         group_means_delta_win.(group{1})=squeeze(mean(delta_win.(group{1}),2));
         group_stes_delta_win.(group{1})=squeeze(std(delta_win.(group{1}),0,2))./sqrt(size(delta_win.(group{1}),2));
         group_means_delta_loss.(group{1})=squeeze(mean(delta_loss.(group{1}),2));
         group_stes_delta_loss.(group{1})=squeeze(std(delta_loss.(group{1}),0,2))./sqrt(size(delta_loss.(group{1}),2));
         
         group_means_diff_delta_win.(group{1})=squeeze(mean(diff_delta_win.(group{1}),2));
         group_stes_diff_delta_win.(group{1})=squeeze(std(diff_delta_win.(group{1}),0,2))./sqrt(size(diff_delta_win.(group{1}),2));
         group_means_diff_delta_loss.(group{1})=squeeze(mean(diff_delta_loss.(group{1}),2));
         group_stes_diff_delta_loss.(group{1})=squeeze(std(diff_delta_loss.(group{1}),0,2))./sqrt(size(diff_delta_loss.(group{1}),2));
         %
         group_means_delta_loss_bal.(group{1})=squeeze(mean(delta_loss_bal.(group{1}),2));
         group_stes_delta_loss_bal.(group{1})=squeeze(std(delta_loss_bal.(group{1}),0,2))./sqrt(size(delta_loss_bal.(group{1}),2));
         group_means_means_pun_delivered_bal.(group{1})=squeeze(mean(means_pun_delivered_bal.(group{1}),2));
         group_stes_means_pun_delivered_bal.(group{1})=squeeze(std(means_pun_delivered_bal.(group{1}),0,2))./sqrt(size(means_pun_delivered_bal.(group{1}),2));
         group_means_means_pun_undelivered_bal.(group{1})=squeeze(mean(means_pun_undelivered_bal.(group{1}),2));
         group_stes_means_pun_undelivered_bal.(group{1})=squeeze(std(means_pun_undelivered_bal.(group{1}),0,2))./sqrt(size(means_pun_undelivered_bal.(group{1}),2));
         group_means_delta_win_bal.(group{1})=squeeze(mean(delta_win_bal.(group{1}),2));
         group_stes_delta_win_bal.(group{1})=squeeze(std(delta_win_bal.(group{1}),0,2))./sqrt(size(delta_win_bal.(group{1}),2));
         group_means_means_rew_delivered_bal.(group{1})=squeeze(mean(means_rew_delivered_bal.(group{1}),2));
         group_stes_means_rew_delivered_bal.(group{1})=squeeze(std(means_rew_delivered_bal.(group{1}),0,2))./sqrt(size(means_rew_delivered_bal.(group{1}),2));
         group_means_means_rew_undelivered_bal.(group{1})=squeeze(mean(means_rew_undelivered_bal.(group{1}),2));
         group_stes_means_rew_undelivered_bal.(group{1})=squeeze(std(means_rew_undelivered_bal.(group{1}),0,2))./sqrt(size(means_rew_undelivered_bal.(group{1}),2));
         group_means_delta_loss_semi.(group{1})=squeeze(mean(delta_loss_semi.(group{1}),2));
         group_stes_delta_loss_semi.(group{1})=squeeze(std(delta_loss_semi.(group{1}),0,2))./sqrt(size(delta_loss_semi.(group{1}),2));
         group_means_means_pun_delivered_semi.(group{1})=squeeze(mean(means_pun_delivered_semi.(group{1}),2));
         group_stes_means_pun_delivered_semi.(group{1})=squeeze(std(means_pun_delivered_semi.(group{1}),0,2))./sqrt(size(means_pun_delivered_semi.(group{1}),2));
         group_means_means_pun_undelivered_semi.(group{1})=squeeze(mean(means_pun_undelivered_semi.(group{1}),2));
         group_stes_means_pun_undelivered_semi.(group{1})=squeeze(std(means_pun_undelivered_semi.(group{1}),0,2))./sqrt(size(means_pun_undelivered_semi.(group{1}),2));
         group_means_delta_win_semi.(group{1})=squeeze(mean(delta_win_semi.(group{1}),2));
         group_stes_delta_win_semi.(group{1})=squeeze(std(delta_win_semi.(group{1}),0,2))./sqrt(size(delta_win_semi.(group{1}),2));
         group_means_means_rew_delivered_semi.(group{1})=squeeze(mean(means_rew_delivered_semi.(group{1}),2));
         group_stes_means_rew_delivered_semi.(group{1})=squeeze(std(means_rew_delivered_semi.(group{1}),0,2))./sqrt(size(means_rew_delivered_semi.(group{1}),2));
         group_means_means_rew_undelivered_semi.(group{1})=squeeze(mean(means_rew_undelivered_semi.(group{1}),2));
         group_stes_means_rew_undelivered_semi.(group{1})=squeeze(std(means_rew_undelivered_semi.(group{1}),0,2))./sqrt(size(means_rew_undelivered_semi.(group{1}),2));
        %diff
         group_means_diff_delta_loss_bal.(group{1})=squeeze(mean(diff_delta_loss_bal.(group{1}),2));
         group_stes_diff_delta_loss_bal.(group{1})=squeeze(std(diff_delta_loss_bal.(group{1}),0,2))./sqrt(size(diff_delta_loss_bal.(group{1}),2));
         group_means_diff_means_pun_delivered_bal.(group{1})=squeeze(mean(diff_means_pun_delivered_bal.(group{1}),2));
         group_stes_diff_means_pun_delivered_bal.(group{1})=squeeze(std(diff_means_pun_delivered_bal.(group{1}),0,2))./sqrt(size(diff_means_pun_delivered_bal.(group{1}),2));
         group_means_diff_means_pun_undelivered_bal.(group{1})=squeeze(mean(diff_means_pun_undelivered_bal.(group{1}),2));
         group_stes_diff_means_pun_undelivered_bal.(group{1})=squeeze(std(diff_means_pun_undelivered_bal.(group{1}),0,2))./sqrt(size(diff_means_pun_undelivered_bal.(group{1}),2));
         group_means_diff_delta_win_bal.(group{1})=squeeze(mean(diff_delta_win_bal.(group{1}),2));
         group_stes_diff_delta_win_bal.(group{1})=squeeze(std(diff_delta_win_bal.(group{1}),0,2))./sqrt(size(diff_delta_win_bal.(group{1}),2));
         group_means_diff_means_rew_delivered_bal.(group{1})=squeeze(mean(diff_means_rew_delivered_bal.(group{1}),2));
         group_stes_diff_means_rew_delivered_bal.(group{1})=squeeze(std(diff_means_rew_delivered_bal.(group{1}),0,2))./sqrt(size(diff_means_rew_delivered_bal.(group{1}),2));
         group_means_diff_means_rew_undelivered_bal.(group{1})=squeeze(mean(diff_means_rew_undelivered_bal.(group{1}),2));
         group_stes_diff_means_rew_undelivered_bal.(group{1})=squeeze(std(diff_means_rew_undelivered_bal.(group{1}),0,2))./sqrt(size(diff_means_rew_undelivered_bal.(group{1}),2));
         group_means_diff_delta_loss_semi.(group{1})=squeeze(mean(diff_delta_loss_semi.(group{1}),2));
         group_stes_diff_delta_loss_semi.(group{1})=squeeze(std(diff_delta_loss_semi.(group{1}),0,2))./sqrt(size(diff_delta_loss_semi.(group{1}),2));
         group_means_diff_means_pun_delivered_semi.(group{1})=squeeze(mean(diff_means_pun_delivered_semi.(group{1}),2));
         group_stes_diff_means_pun_delivered_semi.(group{1})=squeeze(std(diff_means_pun_delivered_semi.(group{1}),0,2))./sqrt(size(diff_means_pun_delivered_semi.(group{1}),2));
         group_means_diff_means_pun_undelivered_semi.(group{1})=squeeze(mean(diff_means_pun_undelivered_semi.(group{1}),2));
         group_stes_diff_means_pun_undelivered_semi.(group{1})=squeeze(std(diff_means_pun_undelivered_semi.(group{1}),0,2))./sqrt(size(diff_means_pun_undelivered_semi.(group{1}),2));
         group_means_diff_delta_win_semi.(group{1})=squeeze(mean(diff_delta_win_semi.(group{1}),2));
         group_stes_diff_delta_win_semi.(group{1})=squeeze(std(diff_delta_win_semi.(group{1}),0,2))./sqrt(size(diff_delta_win_semi.(group{1}),2));
         group_means_diff_means_rew_delivered_semi.(group{1})=squeeze(mean(diff_means_rew_delivered_semi.(group{1}),2));
         group_stes_diff_means_rew_delivered_semi.(group{1})=squeeze(std(diff_means_rew_delivered_semi.(group{1}),0,2))./sqrt(size(diff_means_rew_delivered_semi.(group{1}),2));
         group_means_diff_means_rew_undelivered_semi.(group{1})=squeeze(mean(diff_means_rew_undelivered_semi.(group{1}),2));
         group_stes_diff_means_rew_undelivered_semi.(group{1})=squeeze(std(diff_means_rew_undelivered_semi.(group{1}),0,2))./sqrt(size(diff_means_rew_undelivered_semi.(group{1}),2));

end
%
plot_pupil_vol_vs_stable
plot_pupil_received_each_block
plot_pupil_not_received_each_block
plot_pupil_each_block
barplot_baselines
plot_pupil_received_and_not_for_winv_and_lossv_block
plot_pupil_received_and_not_for_semi_vs_balanced_blocks
plot_pupil_response_diff_vol_sta_to_otherschedule
save('model_free_pupil_results.mat')
%%
aa=mean(means_pun_bslines.reboxetine(2,:,:)-means_pun_bslines.reboxetine(1,:,:),3);
bb=mean(means_pun_bslines.placebo(2,:,:)-means_pun_bslines.placebo(1,:,:),3);
pun_bs_diff=[aa';bb'];
aa=mean(means_rew_bslines.reboxetine(2,:,:)-means_rew_bslines.reboxetine(1,:,:),3);
bb=mean(means_rew_bslines.placebo(2,:,:)-means_rew_bslines.placebo(1,:,:),3);
rew_bs_diff=[aa';bb'];