%logistic regression
%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2019_drug_study/code');
getfolders;

%% load questionnaire data 
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
sublist.placebo=sublist.all(Ques.treatment==2);
% for i=1:length(excnum)
%     for fnames=fieldnames(Ques)'
%         Ques.(fnames{1})(excnum(i))=[];
%     end
% end

tn=80;
blkn=4;
abandontn=10; % abandon first excnum trials when calculating posterior probability
start=[0.5,0.5];

%all matrices for orders of 4 blocks: for 1 for both volatile, 2 for win
%volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];

blkname={'both volatile','win volatile','loss volatile','both stable'};
%% model fitting for all data in data filefolder
for group={'reboxetine','placebo','all'}
    ssublist=sublist.(group{1});
    for v=1:2
            for ss=1:size(ssublist,1)
                    sdatadir=[datadir,num2str(ssublist(ss)),'/'];

                    filename=dir([sdatadir,'*_visit_',num2str(v),'_blktype_*.txt']);

                    data=read_txt([sdatadir,filename.name]);

                    %find block index for each block condition
                    for j=1:blkn
                            blkindex(j)=find(blockorders(data.blktype,:)==j);
                    end

                    for i=1:blkn
                            warning('')
                            information=[data.winpos(((blkindex(i)-1)*tn+1):blkindex(i)*tn),data.losspos(((blkindex(i)-1)*tn+1):blkindex(i)*tn)];
                            choice=data.choice(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            resp=true(tn,1);
                            winchosen=data.winchosen(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            losschosen=data.losschosen(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            outchosen=[winchosen,losschosen];
                            RT=data.RT(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            ifswitch=choice(2:end)-choice(1:end-1);
                            ifswitch=1-abs(ifswitch);
                            oneback=outchosen(1:end-1,:);
                            twoback=[nan(1,2);outchosen(1:end-2,:)];
                            threeback=[nan(2,2);outchosen(1:end-3,:)];
                            fourback=[nan(3,2);outchosen(1:end-4,:)];
                            fiveback=[nan(4,2);outchosen(1:end-5,:)];
                           % sixback=[nan(5,2);outchosen(1:end-6,:)];
                            
                            regmatrix=[oneback,twoback,threeback,fourback,fiveback];
                            opts=statset('glmfit');
                            opts.MaxIter=1000000;
                            betas.(group{1})(v,ss,i,:)=glmfit(regmatrix,ifswitch,'binomial','link','logit');
                            [warnMsg,warnID]=lastwarn;
                            if ~isempty(warnMsg)
                                warnlog.(group{1})(v,ss,blkindex(i))=1;
                            end
                    end
            end
    end
end
%save([datadir,'\beta_estimate.mat'])

%% plot bargarph for volatility effect for each block type for visit 1 for all participants
betas_win_all=squeeze(betas.all(1,:,:,2:2:10));
betas_loss_all=-squeeze(betas.all(1,:,:,3:2:11));
mean_beta_win_all=squeeze(mean(betas_win_all));
mean_beta_loss_all=squeeze(mean(betas_loss_all));
sem_beta_win_all=squeeze(std(betas_win_all))./sqrt(size(betas.all,2));
sem_beta_loss_all=squeeze(std(betas_loss_all))./sqrt(size(betas.all,2));
f1=figure('units','inch','position',[5,2,9,7]);
xtl={'t-1','t-2','t-3','t-4','t-5'};
cond1={'Win','Loss'};
for i=1:blkn
        subplot(2,2,i)
        mean_beta=[mean_beta_win_all(i,:)',mean_beta_loss_all(i,:)'];
        sem_beta=[sem_beta_win_all(i,:)',sem_beta_loss_all(i,:)'];
        H=bar(mean_beta);
        H(1).FaceColor = 'green';
        H(2).FaceColor = 'red';
        set(gca,'XTickLabel',xtl,'FontSize',8);
        ylabel('beta estimates','FontSize',8);
        legend('Win','Loss','AutoUpdate','off');
        title(blkname{i})
        hold on
        pos=1;
        for ii=1:length(cond1)
        for iii=1:length(xtl)
          pos=iii+0.14*(-1)^ii;
          plot([pos,pos],[mean_beta(iii,ii)-sem_beta(iii,ii),mean_beta(iii,ii)+sem_beta(iii,ii)],'-k','LineWidth',1)
        end
        end
        hold on
end
sgtitle('all participants @ visit-1')
saveas(f1,[figdir,'\all_visit_1_logistic_regression_on_switch_comp_blocks.png'])

%% plot bargarph comparing groups for each block*visit*valence
for v=1:2
    for group={'reboxetine','placebo'}
            betas_win.(group{1})=squeeze(betas.(group{1})(v,:,:,2:2:10));
            betas_loss.(group{1})=-squeeze(betas.(group{1})(v,:,:,3:2:11));     
            mean_beta_win.(group{1})=squeeze(mean(betas_win.(group{1})));
            mean_beta_loss.(group{1})=squeeze(mean(betas_loss.(group{1})));
            sem_beta_win.(group{1})=squeeze(std(betas_win.(group{1})))./sqrt(size(betas.(group{1}),2));
            sem_beta_loss.(group{1})=squeeze(std(betas_loss.(group{1})))./sqrt(size(betas.(group{1}),2));
    end

    f2=figure('units','inch','position',[5,2,9,7]);
    xtl={'t-1','t-2','t-3','t-4','t-5'};
    group={'reboxetine','placebo'};
    %win
    for i=1:blkn
            subplot(2,2,i)
            mean_beta=[mean_beta_win.placebo(i,:)',mean_beta_win.reboxetine(i,:)'];
            sem_beta=[sem_beta_win.placebo(i,:)',sem_beta_win.reboxetine(i,:)'];
            H=bar(mean_beta);
            H(1).FaceColor = [0.8 0.9 0.6];
            H(2).FaceColor = 'green';
            set(gca,'XTickLabel',xtl,'FontSize',8);
            ylabel('beta estimates','FontSize',8);
            legend('Placebo','Reboxetine','AutoUpdate','off');
            title(blkname{i})
            hold on
            pos=1;
            for ii=1:length(group)
                    for iii=1:length(xtl)
                      pos=iii+0.14*(-1)^ii;
                      plot([pos,pos],[mean_beta(iii,ii)-sem_beta(iii,ii),mean_beta(iii,ii)+sem_beta(iii,ii)],'-k','LineWidth',1)
                    end
            end
            hold on
    end
    sgtitle(['win outcomes @visit-',num2str(v)])
    saveas(f2,[figdir,'\cmp_groups_visit_',num2str(v),'_win_logistic_regression_on_switch.png'])
        %loss
            f2=figure('units','inch','position',[5,2,9,7]);
    for i=1:blkn
            subplot(2,2,i)
            mean_beta=[mean_beta_loss.placebo(i,:)',mean_beta_loss.reboxetine(i,:)'];
            sem_beta=[sem_beta_loss.placebo(i,:)',sem_beta_loss.reboxetine(i,:)'];
            H=bar(mean_beta);
            H(1).FaceColor = [1 0.6 0.6];
            H(2).FaceColor = 'red';
            set(gca,'XTickLabel',xtl,'FontSize',8);
            ylabel('beta estimates','FontSize',8);
            legend('Placebo','Reboxetine','AutoUpdate','off');
            title(blkname{i})
            hold on
            pos=1;
            for ii=1:length(group)
                    for iii=1:length(xtl)
                      pos=iii+0.14*(-1)^ii;
                      plot([pos,pos],[mean_beta(iii,ii)-sem_beta(iii,ii),mean_beta(iii,ii)+sem_beta(iii,ii)],'-k','LineWidth',1)
                    end
            end
            hold on
    end
    sgtitle(['loss outcomes @visit-',num2str(v)])
    saveas(f2,[figdir,'\cmp_groups_visit_',num2str(v),'_loss_logistic_regression_on_switch.png'])
end
%% plot bargarph comparing groups for visit2 - visit1 each block*valence
 for group={'reboxetine','placebo'}
            betas_win.(group{1})=squeeze(betas.(group{1})(2,:,:,2:2:10)-betas.(group{1})(1,:,:,2:2:10));
            betas_loss.(group{1})=-squeeze(betas.(group{1})(2,:,:,3:2:11)-betas.(group{1})(1,:,:,3:2:11));     
            mean_beta_win.(group{1})=squeeze(mean(betas_win.(group{1})));
            mean_beta_loss.(group{1})=squeeze(mean(betas_loss.(group{1})));
            sem_beta_win.(group{1})=squeeze(std(betas_win.(group{1})))./sqrt(size(betas.(group{1}),2));
            sem_beta_loss.(group{1})=squeeze(std(betas_loss.(group{1})))./sqrt(size(betas.(group{1}),2));
   end

    f2=figure('units','inch','position',[5,2,9,7]);
    xtl={'t-1','t-2','t-3','t-4','t-5'};
    group={'reboxetine','placebo'};
    %win
    for i=1:blkn
            subplot(2,2,i)
            mean_beta=[mean_beta_win.placebo(i,:)',mean_beta_win.reboxetine(i,:)'];
            sem_beta=[sem_beta_win.placebo(i,:)',sem_beta_win.reboxetine(i,:)'];
            H=bar(mean_beta);
            H(1).FaceColor = [0.8 0.9 0.6];
            H(2).FaceColor = 'green';
            set(gca,'XTickLabel',xtl,'FontSize',8);
            ylabel('beta estimates','FontSize',8);
            legend('Placebo','Reboxetine','AutoUpdate','off');
            title(blkname{i})
            hold on
            pos=1;
            for ii=1:length(group)
                    for iii=1:length(xtl)
                      pos=iii+0.14*(-1)^ii;
                      plot([pos,pos],[mean_beta(iii,ii)-sem_beta(iii,ii),mean_beta(iii,ii)+sem_beta(iii,ii)],'-k','LineWidth',1)
                    end
            end
            hold on
    end
    sgtitle('win outcomes @visit2 vs visit1')
    saveas(f2,[figdir,'\cmp_groups_diff_win_logistic_regression.png'])
        %loss
            f2=figure('units','inch','position',[5,2,9,7]);
    for i=1:blkn
            subplot(2,2,i)
            mean_beta=[mean_beta_loss.placebo(i,:)',mean_beta_loss.reboxetine(i,:)'];
            sem_beta=[sem_beta_loss.placebo(i,:)',sem_beta_loss.reboxetine(i,:)'];
            H=bar(mean_beta);
            H(1).FaceColor = [1 0.6 0.6];
            H(2).FaceColor = 'red';
            set(gca,'XTickLabel',xtl,'FontSize',8);
            ylabel('beta estimates','FontSize',8);
            legend('Placebo','Reboxetine','AutoUpdate','off');
            title(blkname{i})
            hold on
            pos=1;
            for ii=1:length(group)
                    for iii=1:length(xtl)
                      pos=iii+0.14*(-1)^ii;
                      plot([pos,pos],[mean_beta(iii,ii)-sem_beta(iii,ii),mean_beta(iii,ii)+sem_beta(iii,ii)],'-k','LineWidth',1)
                    end
            end
            hold on
    end
    sgtitle('loss outcomes @visit2 vs visit1')
    saveas(f2,[figdir,'\cmp_groups_diff_loss_logistic_regression_on_switch.png'])