%% 
clear all;
clc;
%%
cd('D:/2019_drug_study/code');
getfolders;

%% load questionnaire data 
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
 sublist.placebo=sublist.all(Ques.treatment==2);
% %exclude 
% excludesub=[2014,1011];%1021
% for group={'reboxetine','placebo','all'}
%     for i=1:length(excludesub)
%         if isempty(find(sublist.(group{1})==excludesub(i)))
%             continue
%         else
%               excnum(i)=find(sublist.(group{1})==excludesub(i));
%         end
%     sublist.(group{1})(excnum(i))=[];
%     end
% end
%%
addpath("plot_code/")
addpath("model-rescorla_wagner_2lr2b/")
%% experiment information setting
tn=80;
blkn=4;
abandontn=5; % abandon first excnum trials when calculating posterior probability
start=[0.5,0.5];

%all matrices for orders of 4 blocks: for 1 for both volatile, 2 for win
%volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];

blkname={'both volatile','win volatile','loss volatile','both stable'};

%% cal switchprobabilities
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

                    %sort the information and choice made for each block condition
                    choiceside=ones(length(data.choice),1);
                    choiceside(strcmp(data.choiceside,'right'))=0;

                    for i=1:blkn
                            sdata(i).information=[data.winpos(((blkindex(i)-1)*tn+1):blkindex(i)*tn),data.losspos(((blkindex(i)-1)*tn+1):blkindex(i)*tn)];
                            sdata(i).choice=data.choice(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            sdata(i).resp=[false(abandontn,1);true(tn-abandontn,1)];
                            sdata(i).winchosen=data.winchosen(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            sdata(i).losschosen=data.losschosen(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            sdata(i).choiceside=choiceside(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                            sdata(i).RT=data.RT(((blkindex(i)-1)*tn+1):blkindex(i)*tn);
                    end
                    blkindex=[];
                    % save choices in a matrix for later use
                    for i=1:blkn
                            choice_mtx(ss,v,i,:)=sdata(i).choice;
                    end

                    %checking probability of chossing left or shape A in each block
                    for i=1:blkn
                        tmp_pct_leftchoice(i,:)=mean(sdata(i).choiceside(sdata(i).resp));
                        tmp_pct_shape1choice(i,:)=mean(sdata(i).choice(sdata(i).resp));

                        %percentages of winchosen and loss not chosen
                        tmp_pct_winchosen(i,:)=nanmean(sdata(i).winchosen);
                        tmp_pct_losschosen(i,:)=nanmean(sdata(i).losschosen);

                        %calculate swtich probabilty for each block
                        tmp_switchprob(i,:)=cal_switch_prob(sdata(i).information,sdata(i).choice,sdata(i).resp);
                    end
                    pct_leftchoice(ss,:)=tmp_pct_leftchoice;
                    pct_shape1choice(ss,:)=tmp_pct_shape1choice;
                    pct_winchosen.(group{1})(v,ss,:)=tmp_pct_winchosen;
                    pct_lossnotchosen.(group{1})(v,ss,:)=1-tmp_pct_losschosen;
                    tmp_pct_leftchoice=[];
                    tmp_pct_shape1choice=[];
                    tmp_pct_winchosen=[];
                    tmp_pct_losschosen=[];
                    switchprob(ss,:)=tmp_switchprob;
                    clear tmp_switchprob
                    
                   %sort RTs
                   for i=1:blkn
                     rt(ss,i)=cal_RTs(sdata(i).RT,sdata(i).winchosen,sdata(i).losschosen,sdata(i).choice);
                   end
                   
                   money.(group{1})(ss,v)=data.totalmoney(end);
                   blktype.(group{1})(ss,v)=data.blktype;
                   
            end
    %% percentages of win vs loss chosen
    pct_winvsloss_chosen.(group{1})=pct_winchosen.(group{1})-(1-pct_lossnotchosen.(group{1}));
    plot_bargarph_percentages_chosen_win_noloss
    plot_bargarph_percentages_of_choosing_left_shape1
    
    plt_2pct(squeeze(pct_winchosen.(group{1})(v,:,:)),squeeze(pct_lossnotchosen.(group{1})(v,:,:)),...
        'percentages','win chosen','lossnotchosen',group{1},v,'win_loss_perf',figdir)

        %% plot bargarph for switch probabilities
        for i=1:size(switchprob,1)
            for j=1:size(switchprob,2)
                  switchprob_both.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_Both');
                  switchprob_nothing.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_Nothing');
                  switchprob_win.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_win');
                  switchprob_nowin.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_nowin');
                  switchprob_loss.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_loss');
                  switchprob_noloss.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_noloss');
                  switchprob_winalone.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_Winalone');
                  switchprob_lossalone.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_Lossalone');
                  switchprob_same.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_same');
                  switchprob_all.(group{1})(v,i,j)=getfield(switchprob,{i,j},'switchprob_all');
            end
        end
        switchprob_win_psns_tr.(group{1})=asstr(switchprob_nowin.(group{1}))-asstr(switchprob_win.(group{1}));
        switchprob_loss_psns_tr.(group{1})=asstr(switchprob_loss.(group{1}))-asstr(switchprob_noloss.(group{1}));
        switchprob_win_psns.(group{1})=switchprob_nowin.(group{1})-switchprob_win.(group{1});
        switchprob_loss_psns.(group{1})=switchprob_loss.(group{1})-switchprob_noloss.(group{1});
        clear switchprob
       plt_2pct(squeeze(switchprob_win_psns.(group{1})(v,:,:)),squeeze(switchprob_loss_psns.(group{1})(v,:,:)),...
       'Positive stay negative switch','Ppsns-win','Ppsns-loss',group{1},v,'win_loss_psns',figdir)
%         plot_bargarph_switchprob_both_nothing
%         plot_bargarph_switchprob_win_noloss
%         plot_bargarph_switchprob_nowin_loss
%         plot_bargarph_switchprob_PNs
        
        %% get RT values
         for i=1:size(rt,1)
            for j=1:size(rt,2)
                  RT_all.(group{1})(v,i,j)=getfield(rt,{i,j},'all');
                  RT_afterwin.(group{1})(v,i,j)=getfield(rt,{i,j},'afterwin');
                  RT_afternowin.(group{1})(v,i,j)=getfield(rt,{i,j},'afternowin');
                  RT_afterloss.(group{1})(v,i,j)=getfield(rt,{i,j},'afterloss');
                  RT_afternoloss.(group{1})(v,i,j)=getfield(rt,{i,j},'afternoloss');
                  RT_stay.(group{1})(v,i,j)=getfield(rt,{i,j},'stay');
                  RT_switch.(group{1})(v,i,j)=getfield(rt,{i,j},'switch');
            end
         end
        
        end
  
% ifplotscatter=false;
% plot_bargarph_switchprob_dif_PNs
% plot_bargarph_switchprob_dif_positive
% plot_bargarph_switchprob_dif_negative
% plot_bargarph_switchprob_dif_both_nothing
% plot_bargarph_percentages_diff_chosen_win_noloss
% plot_bargarph_percentages_diff_chosen_win_noloss
% plot_bargarph_diff_switchprob_nowin_loss
% plot_bargarph_diff_switchprob_win_noloss
end
%%
plt_barplot_perf_diff_reb_pla
plt_barplot_loss_psns_adpt_othervol_diff_reb_pla
%%
T_pct_winchosen=array2table([squeeze(pct_winchosen.all(1,:,:)),squeeze(pct_winchosen.all(2,:,:))],...
    'VariableNames',{'pct_winchosen_visit1_bothv','pct_winchosen_visit1_winv','pct_winchosen_visit1_lossv','pct_winchosen_visit1_boths',...
    'pct_winchosen_visit2_bothv','pct_winchosen_visit2_winv','pct_winchosen_visit2_lossv','pct_winvchosen_visit2_boths'});
T_pct_winchosen_diff=array2table(squeeze(pct_winchosen.all(2,:,:))-squeeze(pct_winchosen.all(1,:,:)),...
    'VariableNames',{'pct_winchosen_diff_bothv','pct_winchosen_diff_winv','pct_winchosen_diff_lossv','pct_winchosen_diff_boths'});
T_pct_lossnotchosen=array2table([squeeze(pct_lossnotchosen.all(1,:,:)),squeeze(pct_lossnotchosen.all(2,:,:))],...
    'VariableNames',{'pct_lossnotchosen_visit1_bothv','pct_lossnotchosen_visit1_winv','pct_lossnotchosen_visit1_lossv','pct_lossnotchosen_visit1_boths',...
    'pct_lossnotchosen_visit2_bothv','pct_lossnotchosen_visit2_winv','pct_lossnotchosen_visit2_lossv','pct_lossnotvchosen_visit2_boths'});
T_pct_lossnotchosen_diff=array2table(squeeze(pct_lossnotchosen.all(2,:,:))-squeeze(pct_lossnotchosen.all(1,:,:)),...
    'VariableNames',{'pct_lossnotchosen_diff_bothv','pct_lossnotchosen_diff_winv','pct_lossnotchosen_diff_lossv','pct_lossnotchosen_diff_boths'});
T_pct_wvsl=array2table([squeeze(pct_winvsloss_chosen.all(1,:,:)),squeeze(pct_winvsloss_chosen.all(2,:,:))],...
    'VariableNames',{'pct_winvsloss_chosen_visit1_bothv','pct_winvsloss_chosen_visit1_winv','pct_winvsloss_chosen_visit1_lossv','pct_winvsloss_chosen_visit1_boths','pct_winvsloss_chosen_visit2_bothv','pct_winvsloss_chosen_visit2_winv','pct_winvsloss_chosen_visit2_lossv','pct_winvsloss_chosen_visit2_boths'});
T_pct_wvsl_diff=array2table(squeeze(pct_winvsloss_chosen.all(2,:,:)-pct_winvsloss_chosen.all(1,:,:)),...
    'VariableNames',{'pct_winvsloss_chosen_diff_bothv','pct_winvsloss_chosen_diff_winv','pct_winvsloss_chosen_diff_lossv','pct_winvsloss_chosen_diff_boths'});
T_nowin=array2table([squeeze(switchprob_nowin.all(1,:,:)),squeeze(switchprob_nowin.all(2,:,:))],...
    'VariableNames',{'nowin_visit1_bothv','nowin_visit1_winv','nowin_visit1_lossv','nowin_visit1_boths','nowin_visit2_bothv','nowin_visit2_winv','nowin_visit2_lossv','nowin_visit2_boths'});
T_loss=array2table([squeeze(switchprob_loss.all(1,:,:)),squeeze(switchprob_loss.all(2,:,:))],...
    'VariableNames',{'loss_visit1_bothv','loss_visit1_winv','loss_visit1_lossv','loss_visit1_boths','loss_visit2_bothv','loss_visit2_winv','loss_visit2_lossv','loss_visit2_boths'});
T_win=array2table([squeeze(switchprob_win.all(1,:,:)),squeeze(switchprob_win.all(2,:,:))],...
    'VariableNames',{'win_visit1_bothv','win_visit1_winv','win_visit1_lossv','win_visit1_boths','win_visit2_bothv','win_visit2_winv','win_visit2_lossv','win_visit2_boths'});
T_noloss=array2table([squeeze(switchprob_noloss.all(1,:,:)),squeeze(switchprob_noloss.all(2,:,:))],...
    'VariableNames',{'noloss_visit1_bothv','noloss_visit1_winv','noloss_visit1_lossv','noloss_visit1_boths','noloss_visit2_bothv','noloss_visit2_winv','noloss_visit2_lossv','noloss_visit2_boths'});
T_winvsnowin=array2table([squeeze(switchprob_win_psns.all(1,:,:)),squeeze(switchprob_win_psns.all(2,:,:))],...
    'VariableNames',{'win_psns_visit1_bothv','win_psns_visit1_winv','win_psns_visit1_lossv','win_psns_visit1_boths','win_psns_visit2_bothv','win_psns_visit2_winv','win_psns_visit2_lossv','win_psns_visit2_boths'});
T_nolossvsloss=array2table([squeeze(switchprob_loss_psns.all(1,:,:)),squeeze(switchprob_loss_psns.all(2,:,:))],...
    'VariableNames',{'loss_psns_visit1_bothv','loss_psns_visit1_winv','loss_psns_visit1_lossv','loss_psns_visit1_boths','loss_psns_visit2_bothv','loss_psns_visit2_winv','loss_psns_visit2_lossv','loss_psns_visit2_boths'});
T_s=table(Ques.subnum,Ques.treatment,money.all(:,1),money.all(:,2),blktype.all(:,1),blktype.all(:,2),'VariableNames',{'subjectnum','treatment','money_v1','money_v2','blktype_visit1','blktype_visit2'});
T_winvsnowin_diff=array2table(squeeze(switchprob_win_psns.all(2,:,:))-squeeze(switchprob_win_psns.all(1,:,:)),...
    'VariableNames',{'win_psns_diff_bothv','win_psns_diff_winv','win_psns_diff_lossv','win_psns_diff_boths'});
T_nolossvsloss_diff=array2table([squeeze(switchprob_loss_psns.all(2,:,:))-squeeze(switchprob_loss_psns.all(1,:,:))],...
    'VariableNames',{'loss_psns_diff_bothv','loss_psns_diff_winv','loss_psns_diff_lossv','loss_psns_diff_boths'});
T_nowin_diff=array2table(squeeze(switchprob_nowin.all(2,:,:))-squeeze(switchprob_nowin.all(1,:,:)),...
    'VariableNames',{'nowin_diff_bothv','nowin_diff_winv','nowin_diff_lossv','nowin_diff_boths'});
T_loss_diff=array2table(squeeze(switchprob_loss.all(2,:,:))-squeeze(switchprob_loss.all(1,:,:)),...
    'VariableNames',{'loss_diff_bothv','loss_diff_winv','loss_diff_lossv','loss_diff_boths'});
T_win_diff=array2table(squeeze(switchprob_win.all(2,:,:))-squeeze(switchprob_win.all(1,:,:)),...
    'VariableNames',{'win_diff_bothv','win_diff_winv','win_diff_lossv','win_diff_boths'});
T_noloss_diff=array2table(squeeze(switchprob_noloss.all(2,:,:))-squeeze(switchprob_noloss.all(1,:,:)),...
    'VariableNames',{'noloss_diff_bothv','noloss_diff_winv','noloss_diff_lossv','noloss_diff_boths'});
T=[T_s,T_pct_winchosen,T_pct_winchosen_diff,T_pct_lossnotchosen,T_pct_lossnotchosen_diff,T_pct_wvsl,T_pct_wvsl_diff,T_nowin,...
    T_loss,T_win,T_noloss,T_winvsnowin,T_winvsnowin_diff,T_nolossvsloss,T_nolossvsloss_diff,T_nowin_diff,T_loss_diff,T_win_diff,T_noloss_diff];
writetable(T,[datadir,'switchprob_reboxetine'])


%% saving behavior result for spss

T_pct_winchosen_tr=array2table([asstr(squeeze(pct_winchosen.all(1,:,:))),asstr(squeeze(pct_winchosen.all(2,:,:)))],...
    'VariableNames',{'pct_winchosen_tr_visit1_bothv','pct_winchosen_tr_visit1_winv','pct_winchosen_tr_visit1_lossv','pct_winchosen_tr_visit1_boths',...
                              'pct_winchosen_tr_visit2_bothv','pct_winchosen_tr_visit2_winv','pct_winchosen_tr_visit2_lossv','pct_winchosen_tr_visit2_boths'});

T_pct_lossnotchosen_tr=array2table([asstr(squeeze(pct_lossnotchosen.all(1,:,:))),asstr(squeeze(pct_lossnotchosen.all(2,:,:)))],...
    'VariableNames',{'pct_lossnotchosen_tr_visit1_bothv','pct_lossnotchosen_tr_visit1_winv','pct_lossnotchosen_tr_visit1_lossv','pct_lossnotchosen_tr_visit1_boths',...
                               'pct_lossnotchosen_tr_visit2_bothv','pct_lossnotchosen_tr_visit2_winv','pct_lossnotchosen_tr_visit2_lossv','pct_lossnotchosen_tr_visit2_boths'});
T_pct_mean_perf_tr=array2table(transpose([mean(asstr(pct_winchosen.all(1,:,:)),3);mean(asstr(pct_lossnotchosen.all(1,:,:)),3);...
                                mean(asstr(pct_winchosen.all(2,:,:)),3);mean(asstr(pct_lossnotchosen.all(2,:,:)),3);...
                                mean(asstr(pct_winchosen.all(2,:,:)),3)-mean(asstr(pct_winchosen.all(1,:,:)),3);...
                                mean(asstr(pct_lossnotchosen.all(2,:,:)),3)-mean(asstr(pct_lossnotchosen.all(1,:,:)),3)]), ...
                 'VariableNames',{'mean_pct_winchosen_tr_visit1','mean_pct_lossnotchosen_tr_visit1',...
                                'mean_pct_winchosen_tr_visit2','mean_pct_lossnotchosen_tr_visit2',...
                                'mean_pct_winchosen_tr_diff','mean_pct_lossnotchosen_tr_diff'});
T_winvsnowin_tr=array2table([squeeze(switchprob_win_psns_tr.all(1,:,:)),squeeze(switchprob_win_psns_tr.all(2,:,:))],...
    'VariableNames',{'win_psns_tr_visit1_bothv','win_psns_tr_visit1_winv','win_psns_tr_visit1_lossv','win_psns_tr_visit1_boths',...
    'win_psns_tr_visit2_bothv','win_psns_tr_visit2_winv','win_psns_tr_visit2_lossv','win_psns_tr_visit2_boths'});
T_winvsnowin_tr_diff=array2table(squeeze(switchprob_win_psns_tr.all(2,:,:))-squeeze(switchprob_win_psns_tr.all(1,:,:)),...
    'VariableNames',{'win_psns_tr_diff_bothv','win_psns_tr_diff_winv','win_psns_tr_diff_lossv','win_psns_tr_diff_boths'});
T_nolossvsloss_tr=array2table([squeeze(switchprob_loss_psns_tr.all(1,:,:)),squeeze(switchprob_loss_psns_tr.all(2,:,:))],...
    'VariableNames',{'loss_psns_tr_visit1_bothv','loss_psns_tr_visit1_winv','loss_psns_tr_visit1_lossv','loss_psns_tr_visit1_boths','loss_psns_tr_visit2_bothv','loss_psns_tr_visit2_winv','loss_psns_tr_visit2_lossv','loss_psns_tr_visit2_boths'});
T_nolossvsloss_tr_diff=array2table(squeeze(switchprob_loss_psns_tr.all(2,:,:))-squeeze(switchprob_loss_psns_tr.all(1,:,:)),...
    'VariableNames',{'loss_psns_tr_diff_bothv','loss_psns_tr_diff_winv','loss_psns_tr_diff_lossv','loss_psns_tr_diff_boths'});
T_s=table(Ques.subnum,Ques.treatment,money.all(:,1),money.all(:,2),blktype.all(:,1),blktype.all(:,2),'VariableNames',{'subjectnum','treatment','money_v1','money_v2','blktype_visit1','blktype_visit2'});
T=[T_s,T_pct_winchosen_tr,T_pct_lossnotchosen_tr,T_pct_mean_perf_tr,...
    T_winvsnowin_tr,T_winvsnowin_tr_diff,...
    T_nolossvsloss_tr,T_nolossvsloss_tr_diff];
writetable(T,[datadir,'2019_reboxetine_study_behavior_asin_sqrt'])
%%
% i=2;
% aa=rand(25,1)-0.5;
% boxplot(pct_lossnotchosen.placebo(1,:,i),'Notch','on'); hold on;scatter(1+aa./10,transpose(squeeze(pct_lossnotchosen.placebo(1,:,i))),'r','filled')