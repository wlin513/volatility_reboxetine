%% script to analyse the pupil data using OLS
clear all
clc
set(0,'DefaultFigureVisible','on')
% 
cd('D:/2019_drug_study/code/pupil/OLS');
clear all
getfolders
clc
ascdatadir=[datadir,'ascfiles/'];
%load model results
load('R_vol_adpt_4betas_othervol_model.mat');
%%
addpath("D:/2019_drug_study/code/model-rescorla_wagner_2lr1b_bias/")
%%
abandontn=5;%abandon first 5 trails for each block 
nblk=4;
tn=80;
samplerate=500;
start=[0.5,0.5];
%matrices for orders of 4 blocks: for 1 for both volatile, 2 for win volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];
%% load subjects and group information
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
sublist.placebo=sublist.all(Ques.treatment==2);
excludesub=[1019,2013,1011];%2013 a lot of nans in the results;1011 for nans for loss in win vol block
%excluding 2004 for its win suprise corr with win outcome for -0.86/-0.88 ???
for group={'reboxetine','placebo'}
    for i=1:length(excludesub)
        if isempty(find(sublist.(group{1})==excludesub(i)))
            continue
        else
              excnum(i)=find(sublist.(group{1})==excludesub(i));
        end
    sublist.(group{1})(excnum(i))=[];
    R_vol_adpt_othervol_model.winalphas.(group{1})(:,excnum(i),:)=[];
    R_vol_adpt_othervol_model.lossalphas.(group{1})(:,excnum(i),:)=[];
    R_vol_adpt_othervol_model.win_sta_adpt_tr.(group{1})(:,excnum(i),:)=[];
    R_vol_adpt_othervol_model.win_vol_adpt_tr.(group{1})(:,excnum(i),:)=[];
    R_vol_adpt_othervol_model.loss_sta_adpt_tr.(group{1})(:,excnum(i),:)=[];
    R_vol_adpt_othervol_model.loss_vol_adpt_tr.(group{1})(:,excnum(i),:)=[];
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
                file_name=['out_',name,'_visit_',num2str(visit),'.mat'];
                load([ascdatadir,file_name])
                load([ascdatadir,name,'_visit_', num2str(visit),'_pupil_exclude_trials.mat']);
                rew_exclude=rew_exclude'; % to exclude trials where more than 50% of the data were missing
                pun_exclude=pun_exclude';
                % interesting data points
                rew_out=out.rew_outcome;
                loss_out=out.pun_outcome;
    
                %apply high pass filter of 0.1HZ
                [z,p,k] = butter(1,0.1/(samplerate/2),'high');  % first order filter? why?
                [sos,g] = zp2sos(z,p,k);
                Hd = dfilt.df2tsos(sos,g);
                a=reshape(rew_out,[],1);
                b=reshape(loss_out,[],1);
                A=a;A(isnan(a))=0;
                B=b;B(isnan(b))=0;
                AA=filter(Hd,A);
                BB=filter(Hd,B);
                AA(isnan(a))=nan;
                BB(isnan(b))=nan;
                rew_out=reshape(AA,size(rew_out));
                loss_out=reshape(BB,size(loss_out));
    
              %loads the behavior data
                filename=dir([sdatadir,'*_visit_',num2str(visit),'_blktype_*.txt']);
                data=read_txt([sdatadir,filename.name]);
                blkorder=blockorders(data.blktype,:);
                %% baseline subtraction
                rew_bl=nanmean(rew_out(:,501:1000),2); %average the data 1 second preceded the outcome for baseline correction
                loss_bl=nanmean(loss_out(:,501:1000),2);
                REW_bl.(group{1})(visit,sub,:)=rew_bl;
                LOSS_bl.(group{1})(visit,sub,:)=loss_bl;

                rew_out_bl=rew_out-rew_bl;
                loss_out_bl=loss_out-loss_bl;
    
                 both_1st=rew_out;
                 both_1st(data.order==0,:)=loss_out(data.order==0,:); %order=1 win first
                 both_2nd=rew_out;
                 both_2nd(data.order==1,:)=loss_out(data.order==1,:); %order=1 win first

                 BOTH_1st.(group{1})(visit,sub,:,:)=both_1st;
                 BOTH_2nd.(group{1})(visit,sub,:,:)=both_2nd;

                %     %baseline subtraction using the data 1 second preceded the first outcome
                %     both_bl=rew_bl;
                %     both_bl(data.order==0,:)=loss_bl(data.order==0,:); %order=1 win first
                %     rew_out_bl=rew_out-repmat(both_bl,1,size(rew_out,2));
                %     loss_out_bl=loss_out-repmat(both_bl,1,size(loss_out,2));

                % %     %no baseline substraction
                %     rew_out_bl=rew_out;
                %     loss_out_bl=loss_out;  
    
                %% exclude trails  for each trial type
                excludetrials=repmat([true(abandontn,1);false(tn-abandontn,1)],nblk,1);
                rewnantrials=isnan(mean(rew_out_bl,2));
                lossnantrials=isnan(mean(loss_out_bl,2));
                bad_trial_rew_block=rewnantrials|rew_exclude; % not exluding first 5 trials
                bad_trial_loss_block=lossnantrials|pun_exclude;
                bad_trials_rew=excludetrials|rewnantrials|rew_exclude;
                bad_trials_loss=excludetrials|lossnantrials|pun_exclude;

                %% limit the timepoints to actual presented timepoints
                %figure out the dur for each outcome


                dur1=round(abs(data.winresonset-data.lossresonset)*samplerate);
                tmp_Tdiff_out2(:,1)=abs(data.winresonset(1:end-1)-data.fixonset(2:end));
                tmp_Tdiff_out2(:,2)=abs(data.lossresonset(1:end-1)-data.fixonset(2:end));
                Tdiff_out2=min(tmp_Tdiff_out2,[],2);
                Tdiff_out2(80,1)=2;Tdiff_out2(160,1)=2;Tdiff_out2(240,1)=2;Tdiff_out2(320,1)=2;
                dur2=round(Tdiff_out2*samplerate);
                wintimeinclude=false(length(data.choice),size(rew_out,2));
                losstimeinclude=false(length(data.choice),size(loss_out,2));
                orderinc=false(length(data.choice),size(rew_out,2));
                orderinc(data.order==1,:)=true;
                for i=1:size(wintimeinclude,1)
                    dur1inc(i,:)=[true(1,dur1(i)+1000),false(1,size(rew_out,2)-dur1(i)-1000)];
                    dur2inc(i,:)=[true(1,dur2(i)+1000),false(1,size(rew_out,2)-dur2(i)-1000)];
                end

                wintimeinclude(orderinc&dur1inc)=true;%order==1 winfirst
                wintimeinclude((not(orderinc))&dur2inc)=true;
                losstimeinclude((not(orderinc))&dur1inc)=true;%order==0 lossfirst
                losstimeinclude(orderinc&dur1inc)=true;

                rew_out_bl_real=rew_out_bl;
                loss_out_bl_real=loss_out_bl;
                rew_out_bl_real(wintimeinclude==0)=nan;
                loss_out_bl_real(losstimeinclude==0)=nan;
                %% define EVs

                % EV-outcome
                % participant receives a win (1 for received 0 for not-received)
                winout=data.winchosen;
                % participant recieves a loss (1 for received 0 for not-received)
                lossout=data.losschosen;


                %EV-volatility (2 for volatile 1 for stable)
                winvolblk=find(blkorder==1|blkorder==2);
                lossvolblk=find(blkorder==1|blkorder==3);
                winvol=ones(320,1);
                winvol((winvolblk(1)-1)*80+1:winvolblk(1)*80)=2;
                winvol((winvolblk(2)-1)*80+1:winvolblk(2)*80)=2;
                lossvol=ones(320,1);
                lossvol((lossvolblk(1)-1)*80+1:lossvolblk(1)*80)=2;
                lossvol((lossvolblk(2)-1)*80+1:lossvolblk(2)*80)=2;
  
                %EV-ifswitch in the next trial(1 for switch, -1 for stay)
                ifswitch=data.choice(2:end)-data.choice(1:end-1);
                ifswitch(ifswitch~=0)=1;
                ifswitch(ifswitch==0)=-1;
                ifswitch([80,160,240,320])=nan;
    

                %EV-order(1 for current outcome presented first,2 for current outcome presented second)
                winorder=data.order;
                winorder(data.order==0)=2;
                lossorder=data.order+1;
    
                %EV-trial number
                trialnum1=repmat((1:1:80)',4,1);
                trialnum2=(1:320)';
    
                %EV-surprise
                for i=1:4
                    information=[data.winpos(((i-1)*tn+1):i*tn),data.losspos(((i-1)*tn+1):i*tn)];
                    choice=data.choice(((i-1)*tn+1):i*tn);
                    bn=find(blkorder==i);
                    %result_2lr_1b_bias=fit_linked_2lr_beta_add(information,choice, start, abandontn);
                    bel=rescorla_wagner_2lr(information,[R_vol_adpt_othervol_model.winalphas.(group{1})(visit,sub,bn)...
                        R_vol_adpt_othervol_model.lossalphas.(group{1})(visit,sub,bn)],start);

                    winev((i-1)*tn+1:i*tn)=bel(:,1);
                    lossev((i-1)*tn+1:i*tn)=bel(:,2);

                    rewlr.(group{1})(visit,sub,i)=R_vol_adpt_othervol_model.winalphas.(group{1})(visit,sub,bn);
                    losslr.(group{1})(visit,sub,i)=R_vol_adpt_othervol_model.lossalphas.(group{1})(visit,sub,bn);
                end
                winev(data.choice==0)=1-winev(data.choice==0);
                lossev(data.choice==0)=1-lossev(data.choice==0);
                winsurp=data.winchosen-winev';
                losssurp=data.losschosen-lossev';

%                 %
%                 load([datadir,'epsilons.mat'])
%                 win_epsilon=nan(320,1);
%                 loss_epsilon=nan(320,1);
%                 for i=1:4
%                 win_epsilon(((i-1)*tn+1):i*tn)=noisy_model_epsilon1.(group{1})(visit,sub,i);
%                 loss_epsilon(((i-1)*tn+1):i*tn)=noisy_model_epsilon2.(group{1})(visit,sub,i);
%                 end
                
               %% create design matrix for regression   
               othervol_win=losssurp;
               othervol_loss=winsurp;
               othervol_win(winorder==1)=mean(losssurp(winorder==2));
               othervol_loss(lossorder==1)=mean(winsurp(lossorder==2));
               des_mat_rew=[ones(size(rew_out_bl,1),1),winsurp,othervol_win,winvol,winorder,trialnum1];
               des_mat_loss=[ones(size(loss_out_bl,1),1),losssurp,othervol_loss,lossvol,lossorder,trialnum1];

               EVs={'Constant','self-PE','other-PE','self-vol','order','trialnumber block'};

               % create design matrix for regression for blocks  
               des_mat_rew_block=[ones(size(rew_out_bl,1),1),winout,winsurp,winorder,trialnum1];
               des_mat_loss_block=[ones(size(loss_out_bl,1),1),lossout,losssurp,lossorder,trialnum1];

               EVs_block={'Constant','outcome','surprise','order','trial number'};
  
            %% calculate difference traces for outcome and ifswitch
            %difference traces for up to 5s after outcome presented
            windiff.(group{1})(visit,sub,:)=nanmean(rew_out_bl(winout==1,:))-nanmean(rew_out_bl(winout==0,:));
            lossdiff.(group{1})(visit,sub,:)=nanmean(loss_out_bl(lossout==1,:))-nanmean(loss_out_bl(lossout==0,:));

            windiff_switch.(group{1})(visit,sub,:)=nanmean(rew_out_bl(ifswitch==1,:))-nanmean(rew_out_bl(ifswitch==-1,:));
            lossdiff_switch.(group{1})(visit,sub,:)=nanmean(loss_out_bl(ifswitch==1,:))-nanmean(loss_out_bl(ifswitch==-1,:)); 
            % difference traces for timepoints actually presented
            windiff_switch_real.(group{1})(visit,sub,:)=nanmean(rew_out_bl_real(ifswitch==1,:))-nanmean(rew_out_bl_real(ifswitch==-1,:));
            lossdiff_switch_real.(group{1})(visit,sub,:)=nanmean(loss_out_bl_real(ifswitch==1,:))-nanmean(loss_out_bl_real(ifswitch==-1,:)); 

            windiff_real.(group{1})(visit,sub,:)=nanmean(rew_out_bl_real(winout==1,:))-nanmean(rew_out_bl_real(winout==0,:));
            lossdiff_real.(group{1})(visit,sub,:)=nanmean(loss_out_bl_real(lossout==1,:))-nanmean(loss_out_bl_real(lossout==0,:)); 

            bsline_rew_vol_vs_stable.(group{1})(visit,sub,:)=nanmean(rew_out_bl(winvol==2,:))-nanmean(rew_out_bl(winvol==1,:));
            bsline_loss_vol_vs_stable.(group{1})(visit,sub,:)=nanmean(loss_out_bl(lossvol==2,:))-nanmean(loss_out_bl(lossvol==1,:));

            order_rew_2nd_vs_1st.(group{1})(visit,sub,:)=nanmean(rew_out(winorder==2,:))-nanmean(rew_out(winorder==1,:));
            order_loss_2nd_vs_1st.(group{1})(visit,sub,:)=nanmean(loss_out(lossorder==2,:))-nanmean(loss_out(lossorder==1,:));

            % difference traces for each block separately
            for i=1:4
                blockindex=false(320,1);
                bn=find(blkorder==i);
                blockindex(80*(bn-1)+1:80*bn)=true;
                windiff_switch_block_real.(group{1})(visit,sub,i,:)=nanmean(rew_out_bl_real(ifswitch==1 & blockindex,:))-nanmean(rew_out_bl_real(ifswitch==-1 & blockindex,:));
                lossdiff_switch_block_real.(group{1})(visit,sub,i,:)=nanmean(loss_out_bl_real(ifswitch==1 & blockindex,:))-nanmean(loss_out_bl_real(ifswitch==-1 & blockindex,:));

                windiff_switch_block.(group{1})(visit,sub,i,:)=nanmean(rew_out_bl(ifswitch==1 & blockindex,:))-nanmean(rew_out_bl(ifswitch==-1 & blockindex,:));
                lossdiff_switch_block.(group{1})(visit,sub,i,:)=nanmean(loss_out_bl(ifswitch==1 & blockindex,:))-nanmean(loss_out_bl(ifswitch==-1 & blockindex,:));

                windiffblock_real.(group{1})(visit,sub,i,:)=nanmean(rew_out_bl_real(winout==1 & blockindex,:))-nanmean(rew_out_bl_real(winout==0 & blockindex,:));
                lossdiffblock_real.(group{1})(visit,sub,i,:)=nanmean(loss_out_bl_real(lossout==1 & blockindex,:))-nanmean(loss_out_bl_real(lossout==0 & blockindex,:));

                windiffblock.(group{1})(visit,sub,i,:)=nanmean(rew_out_bl(winout==1 & blockindex,:))-nanmean(rew_out_bl(winout==0 & blockindex,:));
                lossdiffblock.(group{1})(visit,sub,i,:)=nanmean(loss_out_bl(lossout==1 & blockindex,:))-nanmean(loss_out_bl(lossout==0 & blockindex,:));

            end

    
%     %plot correlation matrix for each participants

% 
%     f1=figure;
%     tmpRrew=Rrew(sub,:,:);
%     imagesc(squeeze(tmpRrew));
%     colorbar;
% 
%     set(gca,'Xtick',1:size(des_mat_rew,2),'XTickLabel',[ ])
%     set(gca,'Ytick',1:size(des_mat_rew,2),'YTickLabel',[ ])
%     for t=1:size(des_mat_rew,2)
%     text(0,t+1,EVs{t});
%     text(t-0.4,size(des_mat_rew,2)+1,EVs{t});
%     end
%     
%     title(['sub:',num2str(sub),' Win Corr Matrix'])
%     H=findobj(gca,'Type','text');
%     set(H,'Rotation',60); % tilt
%     saveas(f1,[figdir,'win_',num2str(sub),'_designmatrix_othervol.png'])
%     
%     
%     f2=figure;
%     tmpRloss=Rloss(sub,:,:);
%     imagesc(squeeze(tmpRloss));
%     colorbar;
%     
%     set(gca,'Xtick',1:size(des_mat_loss,2),'XTickLabel',[ ])
%     set(gca,'Ytick',1:size(des_mat_loss,2),'YTickLabel',[ ])
%     for t=1:size(des_mat_loss,2)
%     text(0,t+1,EVs{t});
%     text(t-0.4,size(des_mat_rew,2)+1,EVs{t});
%     end
%     title(['sub:',num2str(sub),' Loss Corr Matrix'])
%     H=findobj(gca,'Type','text');
%     set(H,'Rotation',60); % tilt
%     saveas(f1,[figdir,'loss_',num2str(sub),'_designmatrix_othervol.png'])

 %% runing the ols for each block separately (not excluding first 5 trials for each block)
      %specify interactions to consider
       % interaction_block{1}={'outcome','switch'};
        interaction_block=[];
        for j=1:size(interaction_block,2)
            tmp=strcat(interaction_block{j}{1,1},{' * '},interaction_block{j}{1,2});
            EVs_block{end+1}=tmp{1,1};
            interinx_block(j,1)=find(strcmp(string(EVs_block),interaction_block{j}{1,1}));
            interinx_block(j,2)=find(strcmp(string(EVs_block),interaction_block{j}{1,2}));
        end
    for i=1:4
        blockindex=false(320,1);
        bn=find(blkorder==i);
        blockindex(80*(bn-1)+1:80*bn)=true;

        %
        rew_out_bl_block(i).data(:,:)=rew_out_bl((not(bad_trial_rew_block)) & blockindex,:);
        loss_out_bl_block(i).data(:,:)=loss_out_bl((not(bad_trial_loss_block)) & blockindex,:);
        des_mat_rew_4block(i).desmat=des_mat_rew_block((not(bad_trial_rew_block)) & blockindex,:);
        des_mat_loss_4block(i).desmat=des_mat_loss_block((not(bad_trial_loss_block)) & blockindex,:);
        
        %mark switch for the last trial for each block
        ifswitchinx_blk=find(strcmp(string(EVs_block),'switch'));
        if ifswitchinx_blk
            %there is no next trial in the final trial of each block, so mark as mean of the rest 
            des_mat_rew_4block(i).desmat(isnan(des_mat_rew_4block(i).desmat(:,ifswitchinx_blk)),ifswitchinx_blk)=nanmean(des_mat_rew_4block(i).desmat(:,ifswitchinx_blk));
            des_mat_loss_4block(i).desmat(isnan(des_mat_loss_4block(i).desmat(:,ifswitchinx_blk)),ifswitchinx_blk)=nanmean(des_mat_loss_4block(i).desmat(:,ifswitchinx_blk));        
        end
        % normalise all regressors other than the constant
        des_mat_rew_4block(i).desmat(:,2:end)=normalize(des_mat_rew_4block(i).desmat(:,2:end));
        des_mat_loss_4block(i).desmat(:,2:end)=normalize(des_mat_loss_4block(i).desmat(:,2:end));
        
        %calculate interaction regressors
        for j=1:size(interaction_block,2)
            des_mat_rew_4block(i).desmat(:,end+1)=normalize(des_mat_rew_4block(i).desmat(:,interinx_block(j,1)).*des_mat_rew_4block(i).desmat(:,interinx_block(j,2)));
            des_mat_loss_4block(i).desmat(:,end+1)=normalize(des_mat_loss_4block(i).desmat(:,interinx_block(j,1)).*des_mat_loss_4block(i).desmat(:,interinx_block(j,2)));
        end  

        ols_range=[1:3500];%-2s-5s
        % run regression
        [rewbeta_block.(group{1})(visit,sub,i,:,:), rewvar_block.(group{1})(visit,sub,i,:,:)]=ols(rew_out_bl_block(i).data(:,ols_range),des_mat_rew_4block(i).desmat,eye(size(des_mat_rew_4block(i).desmat,2)));
        [lossbeta_block.(group{1})(visit,sub,i,:,:), lossvar_block.(group{1})(visit,sub,i,:,:)]=ols(loss_out_bl_block(i).data(:,ols_range),des_mat_loss_4block(i).desmat,eye(size(des_mat_loss_4block(i).desmat,2)));

    end
    
    clear rew_out_bl_block loss_out_bl_block des_mat_rew_4block des_mat_loss_4block
    %% run ols for all trials
    % get rid of trials with missing data
    bad_trials_rew=bad_trials_rew|lossvol==1;
    bad_trials_loss=bad_trials_loss|winvol==1;
    rew_out_bl(bad_trials_rew,:)=[];
    loss_out_bl(bad_trials_loss,:)=[];
    des_mat_rew(bad_trials_rew,:)=[];
    des_mat_loss(bad_trials_loss,:)=[];
    
    ifswitchinx=find(strcmp(string(EVs),'switch'));
    if ifswitchinx
        %there is no next trial in the final trial of each block, so mark as mean of the rest 
        des_mat_rew(isnan(des_mat_rew(:,ifswitchinx)),ifswitchinx)=nanmean(des_mat_rew(:,ifswitchinx));
        des_mat_loss(isnan(des_mat_loss(:,ifswitchinx)),ifswitchinx)=nanmean(des_mat_loss(:,ifswitchinx));        
    end
    % normalise all regressors other than the constant
    des_mat_rew(:,2:end)=normalize(des_mat_rew(:,2:end));
    des_mat_loss(:,2:end)=normalize(des_mat_loss(:,2:end));

    %specify interactions to consider
    interaction{1}={'self-vol','self-PE'};
    for i=1:size(interaction,2)
        tmp=interaction{i}{1,1};
        for j=2:size(interaction{i},2)
        tmp=strcat(tmp,{' * '},interaction{i}{1,j});
        end
        EVs{length(EVs)+1}=tmp{1,1};
        %index which column in the design matrix to use as interaction variables
        tmpregrew=1;
        tmpregloss=1;
        for j=1:size(interaction{i},2)
        interinx=find(strcmp(string(EVs),interaction{i}{1,j}));  
        tmpregrew=des_mat_rew(:,interinx).*tmpregrew;
        tmpregloss=des_mat_loss(:,interinx).*tmpregloss;
        end
        des_mat_rew(:,end+1)=normalize(tmpregrew);
        des_mat_loss(:,end+1)=normalize(tmpregloss);
    end  
    
    %store corrcoef matrix for each participant
    Rrew.(group{1})(visit,sub,:,:)=corrcoef(des_mat_rew(:,2:end));
    Rloss.(group{1})(visit,sub,:,:)=corrcoef(des_mat_loss(:,2:end));
    
    ols_range=[1:3500];%-2s-5s
    % run regression
    %[rewbeta(sub,:,:) rewvar(sub,:,:)]=ols(rew_out_bl,des_mat_rew,[eye(size(des_mat_rew,2));[ 0 0 0 0 0 1 -1]]);
    %[lossbeta(sub,:,:) lossvar(sub,:,:)]=ols(loss_out_bl,des_mat_loss,[eye(size(des_mat_loss,2));[ 0 0 0 0 0 -1 1]]);
    [rewbeta.(group{1})(visit,sub,:,:), rewvar.(group{1})(visit,sub,:,:)]=ols(rew_out_bl(:,ols_range),des_mat_rew,eye(size(des_mat_rew,2)));
    [lossbeta.(group{1})(visit,sub,:,:), lossvar.(group{1})(visit,sub,:,:)]=ols(loss_out_bl(:,ols_range),des_mat_loss,eye(size(des_mat_loss,2)));

    
    
    
    %%
    trialremain_rew.(group{1})(visit,sub,:)=[str2num(name),80-sum(bad_trials_rew(1:80)),80-sum(bad_trials_rew(81:160)),80-sum(bad_trials_rew(161:240)),80-sum(bad_trials_rew(241:320))];
    trialremain_loss.(group{1})(visit,sub,:)=[str2num(name),80-sum(bad_trials_loss(1:80)),80-sum(bad_trials_loss(81:160)),80-sum(bad_trials_loss(161:240)),80-sum(bad_trials_loss(241:320))];
  
   
       end
    end
end
%% save workspace for future use
save('olsresult_othervol.mat');
%load('olsresult.mat');
%% figures
figdir=[figdir,'ols_results'];
if ~exist( figdir,'dir')
     mkdir(figdir);
end
% %% correlate beta with deppression or anxiety scores

%% plot correlation matrix mean or max

f3=figure('units','inch','position',[0,0,9*2,6]);
subplot(1,2,1)
meanRrew=mean(Rrew.(group{1})(visit,:,:,:),2);
imagesc(squeeze(meanRrew));
colorbar;

set(gca,'Xtick',1:size(des_mat_rew,2),'XTickLabel',[ ])
set(gca,'Ytick',1:size(des_mat_rew,2),'YTickLabel',[ ])

for t=2:size(des_mat_rew,2)
text(0,t,EVs{t});
text(t-1.4,size(des_mat_rew,2),EVs{t});
end

title(['Win mean Corr Matrix'])
H=findobj(gca,'Type','text');
set(H,'Rotation',60); % tilt

subplot(1,2,2)
meanRloss=mean(Rloss.(group{1})(visit,:,:,:),2);
imagesc(squeeze(meanRloss));
colorbar;

set(gca,'Xtick',1:size(des_mat_loss,2),'XTickLabel',[ ])
set(gca,'Ytick',1:size(des_mat_loss,2),'YTickLabel',[ ])
for t=2:size(des_mat_loss,2)
text(0,t,EVs{t});
text(t-1.4,size(des_mat_rew,2),EVs{t});
end
title(['Loss mean Corr Matrix'])
H=findobj(gca,'Type','text');
set(H,'Rotation',60); % tilt
saveas(f3,[figdir,'/mean_designmatrix_othervol.png'])
%%
 colorset=[0,0,0;1,0.75,0.8;1,0.5,0;1,0.84,0;...
     0.74,0.99,0.79;0.03,0.18,0.33;0.12,0.56,1;....
     0.54,0.17,0.89;0.75,0.75,0.75;0.85,0.44,0.84;...
     0.94,0.64,0.38;1,0.49,0.25;0.2,0.3,0.4;0.4,0.5,0.7;...
     0.6,0.2,0.5;0.7,0.3,0.6;0.7,0.3,0.2;0.8,0.3,0.4];
 
 if size(colorset,1)<length(EVs)
     nd=length(EVs)-size(colorset,1);
     for i=1:nd
         colorset=[colorset;rand(1,3)];
     end
 end
% %colorset=rand(length(EVs),3);
% %plot all EVs in one figure
% plotrange=2:3;
% timelab=-2000:2:5000-1;
% 
% f1=figure;
% 
% subplot(2,1,1);
% hold on;
% for i=plotrange
% mean1=squeeze(mean(rewbeta.(group{1})(visit,:,i,:)))';
% sem=squeeze(std(rewbeta.(group{1})(visit,:,i,:)))'./sqrt(size(rewbeta.(group{1})(visit,:,i,:),2));
% sem=sem';
% jbfill(timelab,mean1+sem,mean1-sem,colorset(i,:),colorset(i,:),0,0.1);
% hold on
% plot(timelab,mean1,'color',colorset(i,:));
% hold on
% end
% 
% title('Win');
% %plot(timelab,zeros(size(timelab)),'k--');
% 
% subplot(2,1,2);
% for i=plotrange
% mean1=squeeze(mean(lossbeta.(group{1})(visit,:,i,:)))';
% sem=squeeze(std(lossbeta.(group{1})(visit,:,i,:)))'./sqrt(size(lossbeta.(group{1})(visit,:,i,:),2));
% sem=sem';
% jbfill(timelab,mean1+sem,mean1-sem,colorset(i,:),colorset(i,:),0,0.1);
% hold on
% plot(timelab,mean1,'color',colorset(i,:));
% hold on
% end
% 
% 
% title('Loss');
% %plot(timelab,zeros(size(timelab)),'k--');
% 
% hold on;
% subplot(2,1,1)
% plot(timelab,zeros(size(timelab)),'k--');
% hold on;
% subplot(2,1,2);
% plot(timelab,zeros(size(timelab)),'k--');
% %ylim([-0.25,0.2])'
% legend(EVs(plotrange),'Position',[0.785714285714286 0.419642857142857 0.2 0.117857142857143]);
% saveas(f1,[figdir,group{1},'visit',num2str(visit),'ols_out_all_EVs_othervol.png'])
% 
% 
%% plot each EV
for group={'reboxetine','placebo'}
        for visit=1:2
            for i=2:length(EVs)
            timelab=-1000:2:5000-1;

            f2=figure;
            %for i=plotrange

            subplot(2,1,1);
            hold on;
            plot_ev(rewbeta.(group{1})(visit,:,i,501:end),timelab,['Win:      ',EVs{i}],colorset(i,:))


            subplot(2,1,2);
            plot_ev(lossbeta.(group{1})(visit,:,i,501:end),timelab,['Loss      ',EVs{i}],colorset(i,:))

            hold on;
            subplot(2,1,1)
            plot(timelab,zeros(size(timelab)),'k--');
            hold on;
            subplot(2,1,2);
            plot(timelab,zeros(size(timelab)),'k--');
            %ylim([-0.25,0.2])
            %legend(EVs(i),'Position',[0.785714285714286 0.419642857142857 0.2 0.117857142857143]);
            sgtitle([group{1},': visit ',num2str(visit)])
            saveas(f2,[figdir,'/',strrep(EVs{i},'*','_'),'_',group{1},'_visit_',num2str(visit),'ols_out_EV_','_othervol.png'])
            end
        end
end
%%
%% plot each EV visit2 - visit1
for group={'reboxetine','placebo'}
    for i=2:length(EVs)
    timelab=-1000:2:5000-1;

    f2=figure;
    %for i=plotrange

    subplot(2,1,1);
    plot_ev(rewbeta.(group{1})(2,:,i,501:end)-rewbeta.(group{1})(1,:,i,501:end),timelab,['Win      ',EVs{i}],colorset(i,:))
    hold on;

    subplot(2,1,2);
    plot_ev(lossbeta.(group{1})(2,:,i,501:end)-lossbeta.(group{1})(1,:,i,501:end),timelab,['Loss      ',EVs{i}],colorset(i,:))
    hold on;
    subplot(2,1,1)
    plot(timelab,zeros(size(timelab)),'k--');
    hold on;
    subplot(2,1,2);
    plot(timelab,zeros(size(timelab)),'k--');
    %ylim([-0.25,0.2])
    %legend(EVs(i),'Position',[0.785714285714286 0.419642857142857 0.2 0.117857142857143]);
    sgtitle([group{1},': visit2 vs visit1'])
    saveas(f2,[figdir,'/',strrep(EVs{i},'*','_'),'_',group{1},'_diff_ols_out_EV_','_othervol.png'])
    end
end

% 
% 
% 
% %% plot each block ols result
% for i=2:length(EVs_block)
% timelab=-2000:2:5000-1;
% 
% f4=figure;
% %for i=plotrange
% 
% subplot(2,2,1);
% hold on
% plot_ev_block(rewbeta_block(:,1,i,:),timelab,'Both Volatile','g',2)
% hold on
% plot_ev_block(lossbeta_block(:,1,i,:),timelab,'Both Volatile','r',1)
% 
% subplot(2,2,2);
% hold on
% plot_ev_block(rewbeta_block(:,2,i,:),timelab,'Win Volatile','g',2)
% hold on
% plot_ev_block(lossbeta_block(:,2,i,:),timelab,'Win Volatile','r',1)
% 
% subplot(2,2,3);
% hold on
% plot_ev_block(rewbeta_block(:,3,i,:),timelab,'Loss Volatile','g',2)
% hold on
% plot_ev_block(lossbeta_block(:,3,i,:),timelab,'Loss Volatile','r',1)
% 
% subplot(2,2,4);
% hold on
% plot_ev_block(rewbeta_block(:,4,i,:),timelab,'Both Stable','g',2)
% hold on
% plot_ev_block(lossbeta_block(:,4,i,:),timelab,'Both Stable','r',1)
% 
% 
% hold on;
% subplot(2,2,1)
% plot(timelab,zeros(size(timelab)),'k--');
% hold on;
% subplot(2,2,2);
% plot(timelab,zeros(size(timelab)),'k--');
% hold on;
% subplot(2,2,3)
% plot(timelab,zeros(size(timelab)),'k--');
% hold on;
% subplot(2,2,4);
% plot(timelab,zeros(size(timelab)),'k--');
% 
% legend(['win-',EVs_block{i}],['loss-',EVs_block{i}],'Position',[0.83 0.0 0.1 1])
% 
% saveas(f4,[figdir,'ols_out_EV_block_',strrep(EVs_block{i},'*','_'),'_othervol.png'])
% end
% 
% %% define range to plot
% range=[-2,5];%seconds
% timelab=range(1)*1000:2:range(2)*1000-1;
% samplerange=(range(1)+2)*samplerate+1:(range(2)+2)*samplerate;
% 
% fall=figure;
% 
% subplot(4,1,2)
% meanwindiff_switch=squeeze(mean(windiff_switch));
% sem=std(windiff_switch)/sqrt(size(windiff_switch,1));
% jbfill(timelab,meanwindiff_switch(samplerange)+sem(samplerange),meanwindiff_switch(samplerange)-sem(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff_switch(samplerange),'g')
% hold on
% 
% meanlossdiff_switch=squeeze(mean(lossdiff_switch));
% sem=std(lossdiff_switch)/sqrt(size(lossdiff_switch,1));
% jbfill(timelab,meanlossdiff_switch(samplerange)+sem(samplerange),meanlossdiff_switch(samplerange)-sem(samplerange),'r','r',0,0.1);
% plot(timelab,meanlossdiff_switch(samplerange),'r')
% hold on
% title('Pupil response switch vs stay')
% hold on
% plot(range*1000,[0 0],'k--')
% 
% legend('win-switch','loss-switch')
% 
% subplot(4,1,1)
% meanwindiff=squeeze(mean(windiff));
% sem=std(windiff)/sqrt(size(windiff,1));
% jbfill(timelab,meanwindiff(samplerange)+sem(samplerange),meanwindiff(samplerange)-sem(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff(samplerange),'g')
% hold on
% meanlossdiff=squeeze(mean(lossdiff));
% sem=std(lossdiff)/sqrt(size(lossdiff,1));
% jbfill(timelab,meanlossdiff(samplerange)+sem(samplerange),meanlossdiff(samplerange)-sem(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiff(samplerange),'r')
% title('Pupil response received vs not-received')
% hold on
% plot(range*1000,[0 0],'k--')
% legend('win-received','loss-received')
% %saveas(fall,[figdir,'pupil responses all trials original_othervol.png'])
% subplot(4,1,3)
% meanbsline_rew_vol_vs_stable=squeeze(mean(bsline_rew_vol_vs_stable));
% sem=std(bsline_rew_vol_vs_stable)/sqrt(size(bsline_rew_vol_vs_stable,1));
% jbfill(timelab,meanbsline_rew_vol_vs_stable(samplerange)+sem(samplerange),meanbsline_rew_vol_vs_stable(samplerange)-sem(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanbsline_rew_vol_vs_stable(samplerange),'g')
% hold on
% meanbsline_loss_vol_vs_stable=squeeze(mean(bsline_loss_vol_vs_stable));
% sem=std(bsline_loss_vol_vs_stable)/sqrt(size(bsline_loss_vol_vs_stable,1));
% jbfill(timelab,meanbsline_loss_vol_vs_stable(samplerange)+sem(samplerange),meanbsline_loss_vol_vs_stable(samplerange)-sem(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanbsline_loss_vol_vs_stable(samplerange),'r')
% title('Pupil response vol vs stable')
% hold on
% plot(range*1000,[0 0],'k--')
% legend('win-vol','loss-vol')
% subplot(4,1,4)
% meanorder_rew_2nd_vs_1st=squeeze(mean(order_rew_2nd_vs_1st));
% sem=std(order_rew_2nd_vs_1st)/sqrt(size(order_rew_2nd_vs_1st,1));
% jbfill(timelab,meanorder_rew_2nd_vs_1st(samplerange)+sem(samplerange),meanorder_rew_2nd_vs_1st(samplerange)-sem(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanorder_rew_2nd_vs_1st(samplerange),'g')
% hold on
% meanorder_loss_2nd_vs_1st=squeeze(mean(order_loss_2nd_vs_1st));
% sem=std(order_loss_2nd_vs_1st)/sqrt(size(order_loss_2nd_vs_1st,1));
% jbfill(timelab,meanorder_loss_2nd_vs_1st(samplerange)+sem(samplerange),meanorder_loss_2nd_vs_1st(samplerange)-sem(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanorder_loss_2nd_vs_1st(samplerange),'r')
% title('Pupil response 2nd vs 1st on no-baseline corrected data')
% hold on
% plot(range*1000,[0 0],'k--')
% legend('win-order','loss-order')
% saveas(fall,[figdir,'pupil responses all trials original_othervol.png'])
% 
% %%
% 
% %below hasen't been edited
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% range=[-0.5 2];
% for i=2:length(EVs)
% fbetablock1=figure;
% 
% sem11=squeeze(nanstd(rewbeta_block(:,1,i,:))./sqrt(sum(~isnan(rewbeta_block(:,1,i,:)),1)));
% sem12=squeeze(nanstd(rewbeta_block(:,2,i,:))./sqrt(sum(~isnan(rewbeta_block(:,2,i,:)),1)));
% sem13=squeeze(nanstd(rewbeta_block(:,3,i,:))./sqrt(sum(~isnan(rewbeta_block(:,3,i,:)),1)));
% sem14=squeeze(nanstd(rewbeta_block(:,4,i,:))./sqrt(sum(~isnan(rewbeta_block(:,4,i,:)),1)));
% 
% sem21=squeeze(nanstd(lossbeta_block(:,1,i,:))./sqrt(sum(~isnan(lossbeta_block(:,1,i,:)),1)));
% sem22=squeeze(nanstd(lossbeta_block(:,2,i,:))./sqrt(sum(~isnan(lossbeta_block(:,2,i,:)),1)));
% sem23=squeeze(nanstd(lossbeta_block(:,3,i,:))./sqrt(sum(~isnan(lossbeta_block(:,3,i,:)),1)));
% sem24=squeeze(nanstd(lossbeta_block(:,4,i,:))./sqrt(sum(~isnan(lossbeta_block(:,4,i,:)),1)));
% 
% sem11=sem11'; sem12=sem12';sem13=sem13'; sem14=sem14';
% sem21=sem21'; sem22=sem22';sem23=sem23'; sem24=sem24';
% 
% subplot(2,2,1);
% hold on
% meanrewbeta_block1=squeeze(nanmean(rewbeta_block(:,1,i,:)));
% meanlossbeta_block1=squeeze(nanmean(lossbeta_block(:,1,i,:)));
% meanrewbeta_block1=meanrewbeta_block1';
% meanlossbeta_block1=meanlossbeta_block1';
% jbfill(timelab,meanrewbeta_block1+sem11,meanrewbeta_block1-sem11,'g','g',0,0.1);
% hold on
% plot(timelab,meanrewbeta_block1,'g')
% jbfill(timelab,meanlossbeta_block1+sem21,meanlossbeta_block1-sem21,'r','r',0,0.1);
% hold on
% plot(timelab,meanlossbeta_block1,'r')
% ylabel('parameter estimates')
% title('Both Volatile')
% hold on
% plot(range*1000,[0 0],'k--')
% %xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,2);
% hold on
% meanrewbeta_block2=squeeze(nanmean(rewbeta_block(:,2,i,:)));
% meanlossbeta_block2=squeeze(nanmean(lossbeta_block(:,2,i,:)));
% meanrewbeta_block2=meanrewbeta_block2';
% meanlossbeta_block2=meanlossbeta_block2';
% jbfill(timelab,meanrewbeta_block2+sem12,meanrewbeta_block2-sem12,'g','g',0,0.1);
% hold on
% plot(timelab,meanrewbeta_block2,'g')
% jbfill(timelab,meanlossbeta_block2+sem22,meanlossbeta_block2-sem22,'r','r',0,0.1);
% hold on
% plot(timelab,meanlossbeta_block2,'r')
% ylabel('parameter estimates')
% hold on
% plot(range*1000,[0 0],'k--') 
% title('Win Volatile')
% %xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,3);
% hold on
% meanrewbeta_block3=squeeze(nanmean(rewbeta_block(:,3,i,:)));
% meanlossbeta_block3=squeeze(nanmean(lossbeta_block(:,3,i,:)));
% meanrewbeta_block3=meanrewbeta_block3';
% meanlossbeta_block3=meanlossbeta_block3';
% jbfill(timelab,meanrewbeta_block3+sem13,meanrewbeta_block3-sem13,'g','g',0,0.1);
% hold on
% plot(timelab,meanrewbeta_block3,'g')
% jbfill(timelab,meanlossbeta_block3+sem23,meanlossbeta_block3-sem23,'r','r',0,0.1);
% hold on
% plot(timelab,meanlossbeta_block3,'r')
% ylabel('parameter estimates')
% title('Loss Volatile')
% hold on
% plot(range*1000,[0 0],'k--') 
% %xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,4);
% hold on;
% meanrewbeta_block4=squeeze(nanmean(rewbeta_block(:,4,i,:)));
% meanlossbeta_block4=squeeze(nanmean(lossbeta_block(:,4,i,:)));
% meanrewbeta_block4=meanrewbeta_block4';
% meanlossbeta_block4=meanlossbeta_block4';
% jbfill(timelab,meanrewbeta_block4+sem14,meanrewbeta_block4-sem14,'g','g',0,0.1);
% hold on
% plot(timelab,meanrewbeta_block4,'g')
% jbfill(timelab,meanlossbeta_block4+sem24,meanlossbeta_block4-sem24,'r','r',0,0.1);
% hold on
% plot(timelab,meanlossbeta_block4,'r')
% ylabel('parameter estimates')
% title('Both Stable')
% hold on
% plot(range*1000,[0 0],'k--') 
% %xlim(range*1000)
% legend(['win-',EVs{i}],['loss-',EVs{i}],'Position',[0.83 0.0 0.1 1])
% 
% 
% saveas(fbetablock1,[figdir,'ols-',EVs{i},' each block separately_othervol.png'])
% end
% 
% 
% %%
% f_all_real=figure;
% 
% 
% meanwindiff_switch_real=squeeze(mean(windiff_switch_real));
% meanlossdiff_switch_real=squeeze(mean(lossdiff_switch_real));
% meanwindiff_real=squeeze(mean(windiff_real));
% meanlossdiff_real=squeeze(mean(lossdiff_real));
% 
% sem1=std(windiff_switch_real)./sqrt(sum(~isnan(windiff_switch_real),1));
% sem2=std(lossdiff_switch_real)./sqrt(sum(~isnan(lossdiff_switch_real),1));
% sem3=std(windiff_real)./sqrt(sum(~isnan(windiff_real),1));
% sem4=std(lossdiff_real)./sqrt(sum(~isnan(lossdiff_real),1));
% 
% min1=min([sum(~isnan(sem1)),sum(~isnan(sem2)),sum(~isnan(sem3)),sum(~isnan(sem4))]);
% maxt=floor(min1/samplerate);
% 
% %refine plot range
% range_real=[-2,maxt-2];%seconds
% timelab_real=range_real(1)*1000:2:range_real(2)*1000-1;
% samplerange_real=(range_real(1)+2)*samplerate+1:(range_real(2)+2)*samplerate;
% 
% subplot(2,1,2)
% jbfill(timelab_real,meanwindiff_switch_real(samplerange_real)+sem1(samplerange_real),meanwindiff_switch_real(samplerange_real)-sem1(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_switch_real(samplerange_real),'g')
% hold on
% jbfill(timelab_real,meanlossdiff_switch_real(samplerange_real)+sem2(samplerange_real),meanlossdiff_switch_real(samplerange_real)-sem2(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_switch_real(samplerange_real),'r')
% title('Pupil response switch vs stay')
% hold on
% plot(range_real*1000,[0 0],'k--')
% ylim([-0.2 0.3])
% legend('win-switch','loss-switch')
% 
% subplot(2,1,1)
% jbfill(timelab_real,meanwindiff_real(samplerange_real)+sem3(samplerange_real),meanwindiff_real(samplerange_real)-sem3(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_real(samplerange_real),'g')
% hold on
% jbfill(timelab_real,meanlossdiff_real(samplerange_real)+sem4(samplerange_real),meanlossdiff_real(samplerange_real)-sem4(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_real(samplerange_real),'r')
% title('Pupil response received vs not-received')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% legend('win-received','loss-received')
% saveas(f_all_real,[figdir,'pupil responses all trials real_othervol.png'])
% %%
% f_rec_real=figure;
% 
% sem11=squeeze(std(windiffblock_real(:,1,:))./sqrt(sum(~isnan(windiffblock_real(:,1,:)),1)));
% sem12=squeeze(std(windiffblock_real(:,2,:,:))./sqrt(sum(~isnan(windiffblock_real(:,2,:)),1)));
% sem13=squeeze(std(windiffblock_real(:,3,:))./sqrt(sum(~isnan(windiffblock_real(:,3,:)),1)));
% sem14=squeeze(std(windiffblock_real(:,4,:))./sqrt(sum(~isnan(windiffblock_real(:,4,:)),1)));
% 
% sem21=squeeze(std(lossdiffblock_real(:,1,:))./sqrt(sum(~isnan(lossdiffblock_real(:,1,:)),1)));
% sem22=squeeze(std(lossdiffblock_real(:,2,:))./sqrt(sum(~isnan(lossdiffblock_real(:,2,:)),1)));
% sem23=squeeze(std(lossdiffblock_real(:,3,:))./sqrt(sum(~isnan(lossdiffblock_real(:,3,:)),1)));
% sem24=squeeze(std(lossdiffblock_real(:,4,:))./sqrt(sum(~isnan(lossdiffblock_real(:,4,:)),1)));
% 
% sem11=sem11'; sem12=sem12';sem13=sem13'; sem14=sem14';
% sem21=sem21'; sem22=sem22';sem23=sem23'; sem24=sem24';
% min1=min([sum(~isnan(sem11)),sum(~isnan(sem12)),sum(~isnan(sem13)),sum(~isnan(sem14)) sum(~isnan(sem21)),sum(~isnan(sem22)),sum(~isnan(sem23)),sum(~isnan(sem24))]);
% maxt=floor(min1/samplerate);
% 
% %refine plot range
% range_real=[-2,maxt-2];%seconds
% timelab_real=range_real(1)*1000:2:range_real(2)*1000-1;
% samplerange_real=(range_real(1)+2)*samplerate+1:(range_real(2)+2)*samplerate;
% 
% subplot(2,2,1);
% hold on
% meanwindiffblock1_real=squeeze(nanmean(windiffblock_real(:,1,:)));
% meanlossdiffblock1_real=squeeze(nanmean(lossdiffblock_real(:,1,:)));
% meanwindiffblock1_real=meanwindiffblock1_real';
% meanlossdiffblock1_real=meanlossdiffblock1_real';
% jbfill(timelab_real,meanwindiffblock1_real(samplerange_real)+sem11(samplerange_real),meanwindiffblock1_real(samplerange_real)-sem11(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiffblock1_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiffblock1_real(samplerange_real)+sem21(samplerange_real),meanlossdiffblock1_real(samplerange_real)-sem21(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiffblock1_real(samplerange_real),'r')
% ylabel('pupil response received-not received')
% title('Both Volatile')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% %legend('Win-received','loss-received')
% 
% subplot(2,2,2);
% hold on
% meanwindiffblock2_real=squeeze(nanmean(windiffblock_real(:,2,:)));
% meanlossdiffblock2_real=squeeze(nanmean(lossdiffblock_real(:,2,:)));
% meanwindiffblock2_real=meanwindiffblock2_real';
% meanlossdiffblock2_real=meanlossdiffblock2_real';
% jbfill(timelab_real,meanwindiffblock2_real(samplerange_real)+sem12(samplerange_real),meanwindiffblock2_real(samplerange_real)-sem12(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiffblock2_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiffblock2_real(samplerange_real)+sem22(samplerange_real),meanlossdiffblock2_real(samplerange_real)-sem22(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiffblock2_real(samplerange_real),'r')
% ylabel('pupil response received-not received')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% title('Win Volatile')
% %legend('Win-received','loss-received')
% 
% subplot(2,2,3);
% hold on
% meanwindiffblock3_real=squeeze(nanmean(windiffblock_real(:,3,:)));
% meanlossdiffblock3_real=squeeze(nanmean(lossdiffblock_real(:,3,:)));
% meanwindiffblock3_real=meanwindiffblock3_real';
% meanlossdiffblock3_real=meanlossdiffblock3_real';
% jbfill(timelab_real,meanwindiffblock3_real(samplerange_real)+sem13(samplerange_real),meanwindiffblock3_real(samplerange_real)-sem13(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiffblock3_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiffblock3_real(samplerange_real)+sem23(samplerange_real),meanlossdiffblock3_real(samplerange_real)-sem23(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiffblock3_real(samplerange_real),'r')
% ylabel('pupil response received-not received')
% title('Loss Volatile')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% %legend('Win-received','loss-received')
% 
% subplot(2,2,4);
% hold on;
% meanwindiffblock4_real=squeeze(nanmean(windiffblock_real(:,4,:)));
% meanlossdiffblock4_real=squeeze(nanmean(lossdiffblock_real(:,4,:)));
% meanwindiffblock4_real=meanwindiffblock4_real';
% meanlossdiffblock4_real=meanlossdiffblock4_real';
% jbfill(timelab_real,meanwindiffblock4_real(samplerange_real)+sem14(samplerange_real),meanwindiffblock4_real(samplerange_real)-sem14(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiffblock4_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiffblock4_real(samplerange_real)+sem24(samplerange_real),meanlossdiffblock4_real(samplerange_real)-sem24(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiffblock4_real(samplerange_real),'r')
% ylabel('pupil response received-not received')
% title('Both Stable')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% legend('win-received','loss-received','Position',[0.83 0.0 0.1 1])
% 
% 
% saveas(f_rec_real,[figdir,'pupil response each block separately real_othervol.png'])
% %%
% f11=figure;
% 
% sem11=squeeze(std(windiffblock(:,1,:))./sqrt(sum(~isnan(windiffblock(:,1,:)),1)));
% sem12=squeeze(std(windiffblock(:,2,:))./sqrt(sum(~isnan(windiffblock(:,2,:)),1)));
% sem13=squeeze(std(windiffblock(:,3,:))./sqrt(sum(~isnan(windiffblock(:,3,:)),1)));
% sem14=squeeze(std(windiffblock(:,4,:))./sqrt(sum(~isnan(windiffblock(:,4,:)),1)));
% 
% sem21=squeeze(std(lossdiffblock(:,1,:))./sqrt(sum(~isnan(lossdiffblock(:,1,:)),1)));
% sem22=squeeze(std(lossdiffblock(:,2,:))./sqrt(sum(~isnan(lossdiffblock(:,2,:)),1)));
% sem23=squeeze(std(lossdiffblock(:,3,:))./sqrt(sum(~isnan(lossdiffblock(:,3,:)),1)));
% sem24=squeeze(std(lossdiffblock(:,4,:))./sqrt(sum(~isnan(lossdiffblock(:,4,:)),1)));
% 
% sem11=sem11'; sem12=sem12';sem13=sem13'; sem14=sem14';
% sem21=sem21'; sem22=sem22';sem23=sem23'; sem24=sem24';
% 
% subplot(2,2,1);
% hold on
% meanwindiffblock1=squeeze(nanmean(windiffblock(:,1,:)));
% meanlossdiffblock1=squeeze(nanmean(lossdiffblock(:,1,:)));
% meanwindiffblock1=meanwindiffblock1';
% meanlossdiffblock1=meanlossdiffblock1';
% jbfill(timelab,meanwindiffblock1(samplerange)+sem11(samplerange),meanwindiffblock1(samplerange)-sem11(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiffblock1(samplerange),'g')
% jbfill(timelab,meanlossdiffblock1(samplerange)+sem21(samplerange),meanlossdiffblock1(samplerange)-sem21(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiffblock1(samplerange),'r')
% ylabel('pupil response received-not received')
% title('Both Volatile')
% hold on
% plot(range*1000,[0 0],'k--')
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,2);
% hold on
% meanwindiffblock2=squeeze(nanmean(windiffblock(:,2,:)));
% meanlossdiffblock2=squeeze(nanmean(lossdiffblock(:,2,:)));
% meanwindiffblock2=meanwindiffblock2';
% meanlossdiffblock2=meanlossdiffblock2';
% jbfill(timelab,meanwindiffblock2(samplerange)+sem12(samplerange),meanwindiffblock2(samplerange)-sem12(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiffblock2(samplerange),'g')
% jbfill(timelab,meanlossdiffblock2(samplerange)+sem22(samplerange),meanlossdiffblock2(samplerange)-sem22(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiffblock2(samplerange),'r')
% ylabel('pupil response received-not received')
% hold on
% plot(range*1000,[0 0],'k--') 
% title('Win Volatile')
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,3);
% hold on
% meanwindiffblock3=squeeze(nanmean(windiffblock(:,3,:)));
% meanlossdiffblock3=squeeze(nanmean(lossdiffblock(:,3,:)));
% meanwindiffblock3=meanwindiffblock3';
% meanlossdiffblock3=meanlossdiffblock3';
% jbfill(timelab,meanwindiffblock3(samplerange)+sem13(samplerange),meanwindiffblock3(samplerange)-sem13(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiffblock3(samplerange),'g')
% jbfill(timelab,meanlossdiffblock3(samplerange)+sem23(samplerange),meanlossdiffblock3(samplerange)-sem23(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiffblock3(samplerange),'r')
% ylabel('pupil response received-not received')
% title('Loss Volatile')
% hold on
% plot(range*1000,[0 0],'k--') 
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,4);
% hold on;
% meanwindiffblock4=squeeze(nanmean(windiffblock(:,4,:)));
% meanlossdiffblock4=squeeze(nanmean(lossdiffblock(:,4,:)));
% meanwindiffblock4=meanwindiffblock4';
% meanlossdiffblock4=meanlossdiffblock4';
% jbfill(timelab,meanwindiffblock4(samplerange)+sem14(samplerange),meanwindiffblock4(samplerange)-sem14(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiffblock4(samplerange),'g')
% jbfill(timelab,meanlossdiffblock4(samplerange)+sem24(samplerange),meanlossdiffblock4(samplerange)-sem24(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiffblock4(samplerange),'r')
% ylabel('pupil response received-not received')
% title('Both Stable')
% hold on
% plot(range*1000,[0 0],'k--') 
% xlim(range*1000)
% legend('win-received','loss-received','Position',[0.83 0.0 0.1 1])
% 
% 
% saveas(f11,[figdir,'pupil response each block separately_othervol.png'])
% %%
% fswitchblockreal=figure;
% 
% sem11=squeeze(std(windiff_switch_block_real(:,1,:))./sqrt(sum(~isnan(windiff_switch_block_real(:,1,:)),1)));
% sem12=squeeze(std(windiff_switch_block_real(:,2,:))./sqrt(sum(~isnan(windiff_switch_block_real(:,2,:)),1)));
% sem13=squeeze(std(windiff_switch_block_real(:,3,:))./sqrt(sum(~isnan(windiff_switch_block_real(:,3,:)),1)));
% sem14=squeeze(std(windiff_switch_block_real(:,4,:))./sqrt(sum(~isnan(windiff_switch_block_real(:,4,:)),1)));
% 
% sem21=squeeze(std(lossdiff_switch_block_real(:,1,:))./sqrt(sum(~isnan(lossdiff_switch_block_real(:,1,:)),1)));
% sem22=squeeze(std(lossdiff_switch_block_real(:,2,:))./sqrt(sum(~isnan(lossdiff_switch_block_real(:,2,:)),1)));
% sem23=squeeze(std(lossdiff_switch_block_real(:,3,:))./sqrt(sum(~isnan(lossdiff_switch_block_real(:,3,:)),1)));
% sem24=squeeze(std(lossdiff_switch_block_real(:,4,:))./sqrt(sum(~isnan(lossdiff_switch_block_real(:,4,:)),1)));
% 
% sem11=sem11'; sem12=sem12';sem13=sem13'; sem14=sem14';
% sem21=sem21'; sem22=sem22';sem23=sem23'; sem24=sem24';
% min1=min([sum(~isnan(sem11)),sum(~isnan(sem12)),sum(~isnan(sem13)),sum(~isnan(sem14)) sum(~isnan(sem21)),sum(~isnan(sem22)),sum(~isnan(sem23)),sum(~isnan(sem24))]);
% maxt=floor(min1/samplerate);
% 
% %refine plot range
% range_real=[-0.5,maxt-2];%seconds
% timelab_real=range_real(1)*1000:2:range_real(2)*1000-1;
% samplerange_real=(range_real(1)+2)*samplerate+1:(range_real(2)+2)*samplerate;
% 
% 
% subplot(2,2,1);
% hold on
% meanwindiff_switch_block1_real=squeeze(-nanmean(windiff_switch_block_real(:,1,:)));
% meanlossdiff_switch_block1_real=squeeze(nanmean(lossdiff_switch_block_real(:,1,:)));
% meanwindiff_switch_block1_real=meanwindiff_switch_block1_real';
% meanlossdiff_switch_block1_real=meanlossdiff_switch_block1_real';
% jbfill(timelab_real,meanwindiff_switch_block1_real(samplerange_real)+sem11(samplerange_real),meanwindiff_switch_block1_real(samplerange_real)-sem11(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_switch_block1_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiff_switch_block1_real(samplerange_real)+sem21(samplerange_real),meanlossdiff_switch_block1_real(samplerange_real)-sem21(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_switch_block1_real(samplerange_real),'r')
% ylabel('Pupil response switch vs stay')
% title('Both Volatile')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% 
% subplot(2,2,2);
% hold on
% meanwindiff_switch_block2_real=squeeze(nanmean(windiff_switch_block_real(:,2,:)));
% meanlossdiff_switch_block2_real=squeeze(nanmean(lossdiff_switch_block_real(:,2,:)));
% meanwindiff_switch_block2_real=meanwindiff_switch_block2_real';
% meanlossdiff_switch_block2_real=meanlossdiff_switch_block2_real';
% jbfill(timelab_real,meanwindiff_switch_block2_real(samplerange_real)+sem12(samplerange_real),meanwindiff_switch_block2_real(samplerange_real)-sem12(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_switch_block2_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiff_switch_block2_real(samplerange_real)+sem22(samplerange_real),meanlossdiff_switch_block2_real(samplerange_real)-sem22(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_switch_block2_real(samplerange_real),'r')
% ylabel('Pupil response switch vs stay')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% title('Win Volatile')
% 
% 
% subplot(2,2,3);
% hold on
% meanwindiff_switch_block3_real=squeeze(nanmean(windiff_switch_block_real(:,3,:)));
% meanlossdiff_switch_block3_real=squeeze(nanmean(lossdiff_switch_block_real(:,3,:)));
% meanwindiff_switch_block3_real=meanwindiff_switch_block3_real';
% meanlossdiff_switch_block3_real=meanlossdiff_switch_block3_real';
% jbfill(timelab_real,meanwindiff_switch_block3_real(samplerange_real)+sem13(samplerange_real),meanwindiff_switch_block3_real(samplerange_real)-sem13(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_switch_block3_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiff_switch_block3_real(samplerange_real)+sem23(samplerange_real),meanlossdiff_switch_block3_real(samplerange_real)-sem23(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_switch_block3_real(samplerange_real),'r')
% ylabel('Pupil response switch vs stay')
% title('Loss Volatile')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% 
% 
% subplot(2,2,4);
% hold on;
% meanwindiff_switch_block4_real=squeeze(nanmean(windiff_switch_block_real(:,4,:)));
% meanlossdiff_switch_block4_real=squeeze(nanmean(lossdiff_switch_block_real(:,4,:)));
% meanwindiff_switch_block4_real=meanwindiff_switch_block4_real';
% meanlossdiff_switch_block4_real=meanlossdiff_switch_block4_real';
% jbfill(timelab_real,meanwindiff_switch_block4_real(samplerange_real)+sem14(samplerange_real),meanwindiff_switch_block4_real(samplerange_real)-sem14(samplerange_real),'g','g',0,0.1);
% hold on
% plot(timelab_real,meanwindiff_switch_block4_real(samplerange_real),'g')
% jbfill(timelab_real,meanlossdiff_switch_block4_real(samplerange_real)+sem24(samplerange_real),meanlossdiff_switch_block4_real(samplerange_real)-sem24(samplerange_real),'r','r',0,0.1);
% hold on
% plot(timelab_real,meanlossdiff_switch_block4_real(samplerange_real),'r')
% ylabel('Pupil response switch vs stay')
% title('Both Stable')
% hold on
% plot(range_real*1000,[0 0],'k--') 
% legend('win-switch','loss-switch','Position',[0.83 0.0 0.1 1])
% 
% saveas(fswitchblockreal,[figdir,'pupil response to switch each block separately real_othervol.png'])
% 
% %%
% fswitchblock1=figure;
% 
% sem11=squeeze(std(windiff_switch_block(:,1,:))./sqrt(sum(~isnan(windiff_switch_block(:,1,:)),1)));
% sem12=squeeze(std(windiff_switch_block(:,2,:))./sqrt(sum(~isnan(windiff_switch_block(:,2,:)),1)));
% sem13=squeeze(std(windiff_switch_block(:,3,:))./sqrt(sum(~isnan(windiff_switch_block(:,3,:)),1)));
% sem14=squeeze(std(windiff_switch_block(:,4,:))./sqrt(sum(~isnan(windiff_switch_block(:,4,:)),1)));
% 
% sem21=squeeze(std(lossdiff_switch_block(:,1,:))./sqrt(sum(~isnan(lossdiff_switch_block(:,1,:)),1)));
% sem22=squeeze(std(lossdiff_switch_block(:,2,:))./sqrt(sum(~isnan(lossdiff_switch_block(:,2,:)),1)));
% sem23=squeeze(std(lossdiff_switch_block(:,3,:))./sqrt(sum(~isnan(lossdiff_switch_block(:,3,:)),1)));
% sem24=squeeze(std(lossdiff_switch_block(:,4,:))./sqrt(sum(~isnan(lossdiff_switch_block(:,4,:)),1)));
% 
% sem11=sem11'; sem12=sem12';sem13=sem13'; sem14=sem14';
% sem21=sem21'; sem22=sem22';sem23=sem23'; sem24=sem24';
% 
% subplot(2,2,1);
% hold on
% meanwindiff_switch_block1=squeeze(nanmean(windiff_switch_block(:,1,:)));
% meanlossdiff_switch_block1=squeeze(nanmean(lossdiff_switch_block(:,1,:)));
% meanwindiff_switch_block1=meanwindiff_switch_block1';
% meanlossdiff_switch_block1=meanlossdiff_switch_block1';
% jbfill(timelab,meanwindiff_switch_block1(samplerange)+sem11(samplerange),meanwindiff_switch_block1(samplerange)-sem11(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff_switch_block1(samplerange),'g')
% jbfill(timelab,meanlossdiff_switch_block1(samplerange)+sem21(samplerange),meanlossdiff_switch_block1(samplerange)-sem21(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiff_switch_block1(samplerange),'r')
% ylabel('Pupil response switch vs stay')
% title('Both Volatile')
% hold on
% plot(range*1000,[0 0],'k--')
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,2);
% hold on
% meanwindiff_switch_block2=squeeze(nanmean(windiff_switch_block(:,2,:)));
% meanlossdiff_switch_block2=squeeze(nanmean(lossdiff_switch_block(:,2,:)));
% meanwindiff_switch_block2=meanwindiff_switch_block2';
% meanlossdiff_switch_block2=meanlossdiff_switch_block2';
% jbfill(timelab,meanwindiff_switch_block2(samplerange)+sem12(samplerange),meanwindiff_switch_block2(samplerange)-sem12(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff_switch_block2(samplerange),'g')
% jbfill(timelab,meanlossdiff_switch_block2(samplerange)+sem22(samplerange),meanlossdiff_switch_block2(samplerange)-sem22(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiff_switch_block2(samplerange),'r')
% ylabel('Pupil response switch vs stay')
% hold on
% plot(range*1000,[0 0],'k--') 
% title('Win Volatile')
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,3);
% hold on
% meanwindiff_switch_block3=squeeze(nanmean(windiff_switch_block(:,3,:)));
% meanlossdiff_switch_block3=squeeze(nanmean(lossdiff_switch_block(:,3,:)));
% meanwindiff_switch_block3=meanwindiff_switch_block3';
% meanlossdiff_switch_block3=meanlossdiff_switch_block3';
% jbfill(timelab,meanwindiff_switch_block3(samplerange)+sem13(samplerange),meanwindiff_switch_block3(samplerange)-sem13(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff_switch_block3(samplerange),'g')
% jbfill(timelab,meanlossdiff_switch_block3(samplerange)+sem23(samplerange),meanlossdiff_switch_block3(samplerange)-sem23(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiff_switch_block3(samplerange),'r')
% ylabel('Pupil response switch vs stay')
% title('Loss Volatile')
% hold on
% plot(range*1000,[0 0],'k--') 
% xlim(range*1000)
% %legend('Win-received','loss-received')
% 
% subplot(2,2,4);
% hold on;
% meanwindiff_switch_block4=squeeze(nanmean(windiff_switch_block(:,4,:)));
% meanlossdiff_switch_block4=squeeze(nanmean(lossdiff_switch_block(:,4,:)));
% meanwindiff_switch_block4=meanwindiff_switch_block4';
% meanlossdiff_switch_block4=meanlossdiff_switch_block4';
% jbfill(timelab,meanwindiff_switch_block4(samplerange)+sem14(samplerange),meanwindiff_switch_block4(samplerange)-sem14(samplerange),'g','g',0,0.1);
% hold on
% plot(timelab,meanwindiff_switch_block4(samplerange),'g')
% jbfill(timelab,meanlossdiff_switch_block4(samplerange)+sem24(samplerange),meanlossdiff_switch_block4(samplerange)-sem24(samplerange),'r','r',0,0.1);
% hold on
% plot(timelab,meanlossdiff_switch_block4(samplerange),'r')
% ylabel('Pupil response switch vs stay')
% title('Both Stable')
% hold on
% plot(range*1000,[0 0],'k--') 
% xlim(range*1000)
% legend('win-switch','loss-switch','Position',[0.83 0.0 0.1 1])
% 
% 
% saveas(fswitchblock1,[figdir,'pupil response to switch each block separately_othervol.png'])
% 
% % nf=figure;
% % hold on;
% % plot(timelab,squeeze(mean(windiffblock(:,2,:))),'-b')
% % plot(timelab,squeeze(mean(lossdiffblock(:,2,:))),'-r')
% % plot(timelab,squeeze(mean(windiffblock(:,3,:))),'--b')
% % plot(timelab,squeeze(mean(lossdiffblock(:,3,:))),'--r')
% % legend({'wins win block','losses win block','wins loss block','losses loss block'})
% % 
% % 
% % f2=figure;
% % plot(timelab,squeeze(mean(lossdiffblock(:,2,:)))-squeeze(mean(windiffblock(:,2,:))));
% % hold on;
% % plot(timelab,squeeze(mean(lossdiffblock(:,3,:)))-squeeze(mean(windiffblock(:,3,:))),'r');
% % legend({'Positive Block','Negative_Block'});
% % ylabel('Pupil Dilation to Negative vs Positive Outcome');
% 

