%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2019_drug_study/code');
getfolders;
figdir=[figdir,'combined_control/'];
%% load questionnaire data 
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
sublist.allplacebo=sublist.all(Ques.treatment==2);
% for i=1:length(excnum)
%     for fnames=fieldnames(Ques)'
%         Ques.(fnames{1})(excnum(i))=[];
%     end
% end
%load control participants from rivastigmine study
rdatadir='/home/wlin/Documents/2019_drug_study/data/rivastigmine_control_data/';
% participant info
load([rdatadir,'subjectsGood.mat'])
sublist.allplacebo=[sublist.allplacebo;subjectsGood(subjectsGood(:,2)==4,1)];
sublist.all=[sublist.all;subjectsGood(subjectsGood(:,2)==4,1)];
% %exclude 
% excludesub=[1016,2014,9011];
% for group={'reboxetine','allplacebo','all'}
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
abandontn=10; % abandon first excnum trials when calculating posterior probability
start=[0.5,0.5];

%all matrices for orders of 4 blocks: for 1 for both volatile, 2 for win
%volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];

blkname={'both volatile','win volatile','loss volatile','both stable'};

%% cal switchprobabilities
for group={'reboxetine','allplacebo','all'}
    ssublist=sublist.(group{1});
    for v=1:2
            for ss=1:size(ssublist,1)
                if ssublist(ss)<8000
                    sdatadir=[datadir,num2str(ssublist(ss)),'/'];

                    filename=dir([sdatadir,'*_visit_',num2str(v),'_blktype_*.txt']);

                    data=read_txt([sdatadir,filename.name]);
                else
                     filename=dir([rdatadir,num2str(ssublist(ss)),'_visit_',num2str(v),'_blktype_*.txt']);

                    data=read_txt([rdatadir,filename.name]);
                end
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
                            sdata(i).resp=true(tn,1);
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
            end
            display('both v:')
            sublist.(group{1})(find(pct_shape1choice(:,1)<0.2))
            display('both s:')
            sublist.(group{1})(find(pct_shape1choice(:,4)<0.2))
    %% percentages of win vs loss chosen
    pct_winvsloss_chosen.(group{1})=pct_winchosen.(group{1})-(1-pct_lossnotchosen.(group{1}));
    plot_bargarph_percentages_chosen_win_noloss
    plot_bargarph_percentages_of_choosing_left_shape1
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
        clear switchprob
        plot_bargarph_switchprob_both_nothing
        plot_bargarph_switchprob_win_noloss
        plot_bargarph_switchprob_nowin_loss
        plot_bargarph_switchprob_PNs
        
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
plot_bargarph_switchprob_dif_PNs
plot_bargarph_switchprob_dif_positive
plot_bargarph_switchprob_dif_negative
plot_bargarph_switchprob_dif_both_nothing
plot_bargarph_percentages_diff_chosen_win_noloss
end
%%
group={'allplacebo'};
switchprob_PNs_v1_pla=[squeeze(1-switchprob_win.(group{1})(1,:,:)-switchprob_nowin.(group{1})(1,:,:)),squeeze(1-switchprob_noloss.(group{1})(1,:,:)-switchprob_loss.(group{1})(1,:,:))];
switchprob_PNs_v2_pla=[squeeze(1-switchprob_win.(group{1})(2,:,:)-switchprob_nowin.(group{1})(2,:,:)),squeeze(1-switchprob_noloss.(group{1})(2,:,:)-switchprob_loss.(group{1})(2,:,:))];
switchprob_PNs_v1_pla_win=[squeeze(1-switchprob_win.(group{1})(1,:,:)-switchprob_nowin.(group{1})(1,:,:))];
switchprob_PNs_v2_pla_win=[squeeze(1-switchprob_win.(group{1})(2,:,:)-switchprob_nowin.(group{1})(2,:,:))];
switchprob_PNs_v1_pla_loss=[squeeze(1-switchprob_noloss.(group{1})(1,:,:)-switchprob_loss.(group{1})(1,:,:))];
switchprob_PNs_v2_pla_loss=[squeeze(1-switchprob_noloss.(group{1})(2,:,:)-switchprob_loss.(group{1})(2,:,:))];

switchprob_dif_PN_pla=switchprob_PNs_v1_pla-switchprob_PNs_v2_pla;
switchprob_dif_PN_pla_win=switchprob_PNs_v1_pla_win-switchprob_PNs_v2_pla_win;
switchprob_dif_PN_pla_loss=switchprob_PNs_v1_pla_loss-switchprob_PNs_v2_pla_loss;
group={'reboxetine'};
switchprob_PNs_v1_reb=[squeeze(1-switchprob_win.(group{1})(1,:,:)-switchprob_nowin.(group{1})(1,:,:)),squeeze(1-switchprob_noloss.(group{1})(1,:,:)-switchprob_loss.(group{1})(1,:,:))];
switchprob_PNs_v2_reb=[squeeze(1-switchprob_win.(group{1})(2,:,:)-switchprob_nowin.(group{1})(2,:,:)),squeeze(1-switchprob_noloss.(group{1})(2,:,:)-switchprob_loss.(group{1})(2,:,:))];
switchprob_dif_PN_reb=switchprob_PNs_v1_reb-switchprob_PNs_v2_reb;
switchprob_PNs_v1_reb_win=[squeeze(1-switchprob_win.(group{1})(1,:,:)-switchprob_nowin.(group{1})(1,:,:))];
switchprob_PNs_v2_reb_win=[squeeze(1-switchprob_win.(group{1})(2,:,:)-switchprob_nowin.(group{1})(2,:,:))];
switchprob_PNs_v1_reb_loss=[squeeze(1-switchprob_noloss.(group{1})(1,:,:)-switchprob_loss.(group{1})(1,:,:))];
switchprob_PNs_v2_reb_loss=[squeeze(1-switchprob_noloss.(group{1})(2,:,:)-switchprob_loss.(group{1})(2,:,:))];
switchprob_dif_PN_reb_win=switchprob_PNs_v1_reb_win-switchprob_PNs_v2_reb_win;
switchprob_dif_PN_reb_loss=switchprob_PNs_v1_reb_loss-switchprob_PNs_v2_reb_loss;

[h,p]=ttest(switchprob_dif_PN_pla_win,switchprob_dif_PN_reb_win)
[h,p]=ttest(switchprob_dif_PN_pla_loss,switchprob_dif_PN_reb_loss)


%%
T_pct_wvsl=array2table([squeeze(pct_winvsloss_chosen.all(1,:,:)),squeeze(pct_winvsloss_chosen.all(2,:,:))],...
    'VariableNames',{'pct_winvsloss_chosen_visit1_bothv','pct_winvsloss_chosen_visit1_winv','pct_winvsloss_chosen_visit1_lossv','pct_winvsloss_chosen_visit1_boths','pct_winvsloss_chosen_visit2_bothv','pct_winvsloss_chosen_visit2_winv','pct_winvsloss_chosen_visit2_lossv','pct_winvsloss_chosen_visit2_boths'});
T_pct_wvsl_diff=array2table(squeeze(pct_winvsloss_chosen.all(2,:,:)-pct_winvsloss_chosen.all(1,:,:)),...
    'VariableNames',{'pct_winvsloss_chosen_diff_bothv','pct_winvsloss_chosen_diff_winv','pct_winvsloss_chosen_diff_lossv','pct_winvsloss_chosen_diff_boths'});
switchprob_winPN.all=(1-switchprob_win.all)-switchprob_nowin.all;
switchprob_lossPN.all=(1-switchprob_noloss.all)-switchprob_loss.all;
T_nowin=array2table([squeeze(switchprob_nowin.all(1,:,:)),squeeze(switchprob_nowin.all(2,:,:))],...
    'VariableNames',{'nowin_visit1_bothv','nowin_visit1_winv','nowin_visit1_lossv','nowin_visit1_boths','nowin_visit2_bothv','nowin_visit2_winv','nowin_visit2_lossv','nowin_visit2_boths'});
T_loss=array2table([squeeze(switchprob_loss.all(1,:,:)),squeeze(switchprob_loss.all(2,:,:))],...
    'VariableNames',{'loss_visit1_bothv','loss_visit1_winv','loss_visit1_lossv','loss_visit1_boths','loss_visit2_bothv','loss_visit2_winv','loss_visit2_lossv','loss_visit2_boths'});
T_win=array2table([squeeze(switchprob_win.all(1,:,:)),squeeze(switchprob_win.all(2,:,:))],...
    'VariableNames',{'win_visit1_bothv','win_visit1_winv','win_visit1_lossv','win_visit1_boths','win_visit2_bothv','win_visit2_winv','win_visit2_lossv','win_visit2_boths'});
T_noloss=array2table([squeeze(switchprob_noloss.all(1,:,:)),squeeze(switchprob_noloss.all(2,:,:))],...
    'VariableNames',{'noloss_visit1_bothv','noloss_visit1_winv','noloss_visit1_lossv','noloss_visit1_boths','noloss_visit2_bothv','noloss_visit2_winv','noloss_visit2_lossv','noloss_visit2_boths'});
T_winvsnowin=array2table([squeeze(switchprob_winPN.all(1,:,:)),squeeze(switchprob_winPN.all(2,:,:))],...
    'VariableNames',{'winPN_visit1_bothv','winPN_visit1_winv','winPN_visit1_lossv','winPN_visit1_boths','winPN_visit2_bothv','winPN_visit2_winv','winPN_visit2_lossv','winPN_visit2_boths'});
T_nolossvsloss=array2table([squeeze(switchprob_lossPN.all(1,:,:)),squeeze(switchprob_lossPN.all(2,:,:))],...
    'VariableNames',{'lossPN_visit1_bothv','lossPN_visit1_winv','lossPN_visit1_lossv','lossPN_visit1_boths','lossPN_visit2_bothv','lossPN_visit2_winv','lossPN_visit2_lossv','lossPN_visit2_boths'});
T_s=table(Ques.subnum,Ques.treatment,money.all(:,1),money.all(:,2),'VariableNames',{'subjectnum','treatment','money_v1','money_v2'});
T=[T_s,T_pct_wvsl,T_pct_wvsl_diff,T_nowin,T_loss,T_win,T_noloss,T_winvsnowin,T_nolossvsloss];
writetable(T,[datadir,'switchprob_reboxetine'])
%%
as_pct_wvsl_v1=asin(sqrt(pct_winchosen.all(1,:,:)))-asin(sqrt(1-pct_lossnotchosen.all(1,:,:)));
as_pct_wvsl_v2=asin(sqrt(pct_winchosen.all(2,:,:)))-asin(sqrt(1-pct_lossnotchosen.all(2,:,:)));
T_pct_wvsl=array2table([squeeze(as_pct_wvsl_v1),squeeze(as_pct_wvsl_v2)],...
    'VariableNames',{'pct_winvsloss_visit1_bothv','pct_winvsloss_chosen_visit1_winv','pct_winvsloss_chosen_visit1_lossv','pct_winvsloss_chosen_visit1_boths','pct_winvsloss_chosen_visit2_bothv','pct_winvsloss_chosen_visit2_winv','pct_winvsloss_chosen_visit2_lossv','pct_winvsloss_chosen_visit2_boths'});
T_pct_wvsl_diff=array2table(squeeze(as_pct_wvsl_v2-as_pct_wvsl_v1),...
    'VariableNames',{'pct_winvsloss_chosen_diff_bothv','pct_winvsloss_chosen_diff_winv','pct_winvsloss_chosen_diff_lossv','pct_winvsloss_chosen_diff_boths'});
%
T_nowin=array2table([asin(sqrt(squeeze(switchprob_nowin.all(1,:,:)))),asin(sqrt(squeeze(switchprob_nowin.all(2,:,:))))],...
    'VariableNames',{'nowin_visit1_bothv','nowin_visit1_winv','nowin_visit1_lossv','nowin_visit1_boths','nowin_visit2_bothv','nowin_visit2_winv','nowin_visit2_lossv','nowin_visit2_boths'});
T_nowin_diff=array2table(asin(sqrt(squeeze(switchprob_nowin.all(2,:,:))))-asin(sqrt(squeeze(switchprob_nowin.all(1,:,:)))),...
    'VariableNames',{'nowin_diff_bothv','nowin_diff_winv','nowin_diff_lossv','nowin_diff_boths'});
%
T_loss=array2table([asin(sqrt(squeeze(switchprob_loss.all(1,:,:)))),asin(sqrt(squeeze(switchprob_loss.all(2,:,:))))],...
    'VariableNames',{'loss_visit1_bothv','loss_visit1_winv','loss_visit1_lossv','loss_visit1_boths','loss_visit2_bothv','loss_visit2_winv','loss_visit2_lossv','loss_visit2_boths'});
T_loss_diff=array2table(asin(sqrt(squeeze(switchprob_loss.all(2,:,:))))-asin(sqrt(squeeze(switchprob_loss.all(1,:,:)))),...
    'VariableNames',{'loss_diff_bothv','loss_diff_winv','loss_diff_lossv','loss_diff_boths'});
%
T_win=array2table([asin(sqrt(squeeze(switchprob_win.all(1,:,:)))),asin(sqrt(squeeze(switchprob_win.all(2,:,:))))],...
    'VariableNames',{'win_visit1_bothv','win_visit1_winv','win_visit1_lossv','win_visit1_boths','win_visit2_bothv','win_visit2_winv','win_visit2_lossv','win_visit2_boths'});
T_win_diff=array2table(asin(sqrt(squeeze(switchprob_win.all(2,:,:))))-asin(sqrt(squeeze(switchprob_win.all(1,:,:)))),...
    'VariableNames',{'win_diff_bothv','win_diff_winv','win_diff_lossv','win_diff_boths'});
%
T_noloss=array2table([asin(sqrt(squeeze(switchprob_noloss.all(1,:,:)))),asin(sqrt(squeeze(switchprob_noloss.all(2,:,:))))],...
    'VariableNames',{'noloss_visit1_bothv','noloss_visit1_winv','noloss_visit1_lossv','noloss_visit1_boths','noloss_visit2_bothv','noloss_visit2_winv','noloss_visit2_lossv','noloss_visit2_boths'});
T_noloss_diff=array2table(asin(sqrt(squeeze(switchprob_noloss.all(2,:,:))))-asin(sqrt(squeeze(switchprob_noloss.all(1,:,:)))),...
    'VariableNames',{'noloss_diff_bothv','noloss_diff_winv','noloss_diff_lossv','noloss_diff_boths'});
%
switchprob_winPN.all=asin(sqrt((1-switchprob_win.all)))-asin(sqrt(switchprob_nowin.all));
switchprob_lossPN.all=asin(sqrt((1-switchprob_noloss.all)))-asin(sqrt(switchprob_loss.all));
T_winvsnowin=array2table([squeeze(switchprob_winPN.all(1,:,:)),squeeze(switchprob_winPN.all(2,:,:))],...
    'VariableNames',{'winPN_visit1_bothv','winPN_visit1_winv','winPN_visit1_lossv','winPN_visit1_boths','winPN_visit2_bothv','winPN_visit2_winv','winPN_visit2_lossv','winPN_visit2_boths'});
T_winvsnowin_diff=array2table(squeeze(switchprob_winPN.all(2,:,:))-squeeze(switchprob_winPN.all(1,:,:)),...
    'VariableNames',{'winPN_diff_bothv','winPN_diff_winv','winPN_diff_lossv','winPN_diff_boths'});
T_nolossvsloss=array2table([squeeze(switchprob_lossPN.all(1,:,:)),squeeze(switchprob_lossPN.all(2,:,:))],...
    'VariableNames',{'lossPN_visit1_bothv','lossPN_visit1_winv','lossPN_visit1_lossv','lossPN_visit1_boths','lossPN_visit2_bothv','lossPN_visit2_winv','lossPN_visit2_lossv','lossPN_visit2_boths'});
T_nolossvsloss_diff=array2table(squeeze(switchprob_lossPN.all(2,:,:))-squeeze(switchprob_lossPN.all(1,:,:)),...
    'VariableNames',{'lossPN_diff_bothv','lossPN_diff_winv','lossPN_diff_lossv','lossPN_diff_boths'});
T_s=table(Ques.subnum,Ques.treatment,money.all(:,1),money.all(:,2),'VariableNames',{'subjectnum','treatment','money_v1','money_v2'});
T=[T_s,T_pct_wvsl,T_pct_wvsl_diff,T_nowin,T_nowin_diff,T_loss,T_loss_diff,T_win,T_win_diff,T_noloss,T_noloss_diff,T_winvsnowin,T_winvsnowin_diff,T_nolossvsloss,T_nolossvsloss_diff];
writetable(T,[datadir,'switchprob_reboxetine_asin_sqrt'])

%% stat tests for RT

T_rt_all=array2table([squeeze(RT_all.all(1,:,:)),squeeze(RT_all.all(2,:,:))],...
    'VariableNames',{'v1_bothv','v1_winv','v1_lossv','v1_boths','v2_bothv','v2_winv','v2_lossv','v2_boths'});
T_RT_all=[T_s,T_rt_all];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
     'VariableNames',{'visit','winvol','lossvol'});
rm=fitrm(T_RT_all,'v1_bothv-v2_boths ~ treatment','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','visit+winvol+lossvol')
writetable(T_RT_all,[datadir,'reactiontime'])

T_rt_all=array2table([squeeze(RT_all.all(1,:,:)),squeeze(RT_all.all(2,:,:))],...
    'VariableNames',{'v1_bothv','v1_winv','v1_lossv','v1_boths','v2_bothv','v2_winv','v2_lossv','v2_boths'});
T_RT_all=[T_s,T_rt_all];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'1';'2';'2';'3'},2,1),...
     'VariableNames',{'visit','overvol'});
rm=fitrm(T_RT_all,'v1_bothv-v2_boths ~ treatment','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','visit+overvol')
%%
T_rt_all=array2table([squeeze(RT_switch.all(1,:,:)),squeeze(RT_switch.all(2,:,:))],...
    'VariableNames',{'v1_bothv','v1_winv','v1_lossv','v1_boths','v2_bothv','v2_winv','v2_lossv','v2_boths'});
T_RT_all=[T_s,T_rt_all];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'1';'2';'2';'3'},2,1),...
     'VariableNames',{'visit','overvol'});
rm=fitrm(T_RT_all,'v1_bothv-v2_boths ~ treatment','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','visit+overvol')
writetable(T_RT_all,[datadir,'reactiontime_switch'])
%%
T_pct_winvsloss=[T_s,T_pct_wvsl];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
     'VariableNames',{'visit','winvol','lossvol'});
 rm=fitrm(T_pct_winvsloss,'pct_winvsloss_visit1_bothv-pct_winvsloss_chosen_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl=ranova(rm,'WithinModel','visit*winvol*lossvol')
 %%
 T_pct_neg=[T_s,T_nowin,T_loss];
within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
     'VariableNames',{'valence','visit','winvol','lossvol'});
 rm=fitrm(T_pct_neg,'nowin_visit1_bothv-loss_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')

 %%
  T_pct_neg=[T_s,T_win,T_noloss];
within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
     'VariableNames',{'valence','visit','winvol','lossvol'});
 rm=fitrm(T_pct_neg,'win_visit1_bothv-noloss_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
%%
T_pct_pn=[T_s,T_winvsnowin,T_nolossvsloss];
within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
     'VariableNames',{'valence','visit','winvol','lossvol'});
 rm=fitrm(T_pct_pn,'winPN_visit1_bothv-lossPN_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
%%
for i=1:4
scatter(i*ones(1,length(sublist.allplacebo)),squeeze(switchprob_nowin.allplacebo(2,:,i)-switchprob_nowin.allplacebo(1,:,i)))
hold on
end
sublist.all(find(switchprob_PNs(:,8)>0.6))