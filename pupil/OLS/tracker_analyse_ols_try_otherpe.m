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
               des_mat_rew=[ones(size(rew_out_bl,1),1),winsurp,othervol_win,winvol,lossvol,winorder,trialnum1];
               des_mat_loss=[ones(size(loss_out_bl,1),1),losssurp,othervol_loss,lossvol,winvol,lossorder,trialnum1];

               EVs={'Constant','self-PE','other-PE','self-vol','other-vol','order','trialnumber block'};

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
%     saveas(f1,[figdir,'win_',num2str(sub),'_designmatrix.png'])
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
%     saveas(f1,[figdir,'loss_',num2str(sub),'_designmatrix.png'])

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
    interaction{2}={'other-vol','self-PE'};
    interaction{3}={'self-vol','other-vol','self-PE'};
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
save('olsresult.mat');
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
saveas(f3,[figdir,'/mean_designmatrix.png'])
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
% saveas(f1,[figdir,group{1},'visit',num2str(visit),'ols_out_all_EVs.png'])
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
            plot_ev(rewbeta.(group{1})(visit,:,i,501:end),timelab,['Win:      ',EVs{i}],'',colorset(i,:))


            subplot(2,1,2);
            plot_ev(lossbeta.(group{1})(visit,:,i,501:end),timelab,['Loss:      ',EVs{i}],'',colorset(i,:))

            hold on;
            subplot(2,1,1)
            plot(timelab,zeros(size(timelab)),'k--');
            hold on;
            subplot(2,1,2);
            plot(timelab,zeros(size(timelab)),'k--');
            %ylim([-0.25,0.2])
            %legend(EVs(i),'Position',[0.785714285714286 0.419642857142857 0.2 0.117857142857143]);
            sgtitle([group{1},': visit ',num2str(visit)])
            saveas(f2,[figdir,'/',strrep(EVs{i},'*','_'),'_',group{1},'_visit_',num2str(visit),'ols_out_EV_','.png'])
            end
        end
end
%%
plt_pupil_loss_selfv_otherv_pe
plt_pupil_rew_selfv_otherv_pe
plt_pupil_rew_pe
plt_pupil_loss_pe
%%
f=figure;
plot_ev(lossbeta.reboxetine(1,:,2,501:end),timelab,'','Parameter estimates for loss PE',[0.3 0.9, 1],true,1)
hold on
plot_ev(lossbeta.reboxetine(2,:,2,501:end),timelab,'Reboxetine','Parameter estimates for loss PE',[0 0.5, 1],true,2)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_PE_reboxetine.png'])

f=figure;
plot_ev(lossbeta.placebo(1,:,2,501:end),timelab,'','Parameter estimates for loss PE',[1,0.95, 0.020],true,1)
hold on
plot_ev(lossbeta.placebo(2,:,2,501:end),timelab,'Placebo','Parameter estimates for loss PE',[0.937,0.702, 0.322],true,2)
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('visit 1','visit 2')
saveas(f,[figdir,'loss_PE_placebo.png'])

f=figure;
plot_ev(lossbeta.reboxetine(2,:,2,501:end)-lossbeta.reboxetine(1,:,2,501:end),timelab,'','Parameter estimates visit2 vs. visit1',[0.188 0.663, 0.871])
hold on
plot_ev(lossbeta.placebo(2,:,2,501:end)-lossbeta.placebo(1,:,2,501:end),timelab,'','Parameter estimates visit2 vs. visit1',[0.937,0.863, 0.020])
hold on
plot(timelab,zeros(size(timelab)),'k--');
legend('reboxetine','placebo')
saveas(f,[figdir,'loss_PE_reboxetine_v2minusv1_vs_placebo.png'])
% a=mean(lossbeta.reboxetine(2,:,2,1001:end)-lossbeta.reboxetine(1,:,2,1001:end),4);
% b=mean(lossbeta.placebo(2,:,2,1001:end)-lossbeta.placebo(1,:,2,1001:end),4);
% [h,p]=ttest2(a,b)
% a=mean(rewbeta.reboxetine(2,:,2,1001:end)-rewbeta.reboxetine(1,:,2,1001:end),4);
% b=mean(rewbeta.placebo(2,:,2,1001:end)-rewbeta.placebo(1,:,2,1001:end),4);
% [h,p]=ttest2(a,b)
%% plot each EV visit2 - visit1
for group={'reboxetine','placebo'}
    for i=2:length(EVs)
    timelab=-1000:2:5000-1;

    f2=figure;
    %for i=plotrange

    subplot(2,1,1);
    plot_ev(rewbeta.(group{1})(2,:,i,501:end)-rewbeta.(group{1})(1,:,i,501:end),timelab,['Win:      ',EVs{i}],colorset(i,:))
    hold on;

    subplot(2,1,2);
    plot_ev(lossbeta.(group{1})(2,:,i,501:end)-lossbeta.(group{1})(1,:,i,501:end),timelab,['Loss:      ',EVs{i}],colorset(i,:))
    hold on;
    subplot(2,1,1)
    plot(timelab,zeros(size(timelab)),'k--');
    hold on;
    subplot(2,1,2);
    plot(timelab,zeros(size(timelab)),'k--');
    %ylim([-0.25,0.2])
    %legend(EVs(i),'Position',[0.785714285714286 0.419642857142857 0.2 0.117857142857143]);
    sgtitle([group{1},': visit2 vs visit1'])
    saveas(f2,[figdir,'/',strrep(EVs{i},'*','_'),'_',group{1},'_diff_ols_out_EV_','.png'])
    end
end
% 
