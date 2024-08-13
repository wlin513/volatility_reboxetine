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

%%
addpath("plot_code/")
addpath("model-rescorla_wagner_1lr1b/")
addpath("model-rescorla_wagner_2lr1b/")
addpath("model-rescorla_wagner_2lr1b_choice_stickness/")
addpath("model-rescorla_wagner_2lr2b/")
addpath("model-rescorla_wagner_PNPE_2lr1b/")
addpath("model-rescorla_wagner_2lr1b_bias/")
addpath("model-rescorla_wagner_PNPE_2lr1b_opt1/")    
addpath("model-rescorla_wagner_PNPE_2lr2b/")
addpath("model-rescorla_wagner_PNPE_2lr2b_chosen/")
addpath("model-rescorla_wagner_4lr1b_chosen/")
%% experiment information setting
tn=80;
blkn=4;
abandontn=10; % abandon first excnum trials when calculating posterior probability
start=[0.5,0.5];

%matrices for orders of 4 blocks: for 1 for both volatile, 2 for win volatile/loss stable, 3 for loss volatile/win stable, 4 for both stable
blockorders=[1,2,3,4; 1,3,2,4;4,2,3,1;4,3,2,1;2,1,4,3;2,4,1,3;3,4,1,2;3,1,4,2];

blkname={'both volatile','win volatile','loss volatile','both stable'};
%% vol_adpt model from R
%read data
rstandir2='D:/2019_drug_study/data/rstan/RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_4betas_allblocks_onefit/';
file=dir([rstandir2,'*.mat']);
R_vol_adpt_othervol_modeltmp=load([rstandir2,file.name]);
blkname2={'both volatile','win volatile','loss volatile','both stable'};
variablenames={'winalphas','lossalphas','betas'};
for i=1:length(variablenames)
    j=1;
  for blk=blkname2
    R_vol_adpt_othervol_modeltmp2.(variablenames{i})(:,j)=R_vol_adpt_othervol_modeltmp.(variablenames{i}).values(strcmp(blk{1},R_vol_adpt_othervol_modeltmp.(variablenames{i}).block));
    j=j+1;
  end
end
%
variablenames={'win_vol_adpt_tr','loss_vol_adpt_tr','win_sta_adpt_tr','loss_sta_adpt_tr'};
for i=1:length(variablenames)
    R_vol_adpt_othervol_modeltmp2.(variablenames{i})=R_vol_adpt_othervol_modeltmp.(variablenames{i}).values(strcmp('na',R_vol_adpt_othervol_modeltmp.(variablenames{i}).block));
end
%follow the structure of "sublist" leaving excluded data as nan
variablenames={'winalphas','lossalphas','betas','win_vol_adpt_tr','loss_vol_adpt_tr','win_sta_adpt_tr','loss_sta_adpt_tr'};
incsubs=R_vol_adpt_othervol_modeltmp.win_vol_adpt_tr.IDs;
for i=1:length(variablenames)
    for j=1:length(incsubs)
        for group={'all','reboxetine','placebo'}
            subidx=find(sublist.(group{1})==str2double(extractBefore(incsubs{j},'_')));
            if subidx
                visitn=str2double(extractAfter(incsubs{j},'visit_'));
                R_vol_adpt_othervol_model.(variablenames{i}).(group{1})(visitn,subidx,:)=R_vol_adpt_othervol_modeltmp2.(variablenames{i})(j,:);
            end
            clear subidx
        end
    end
end
% %exclude 1007 for visit2 data have R>1.1, exclude 2004 for not being able to converge
% excludesubs=[1029];
% for i=1:length(excludesubs)
%     for group={'all','reboxetine','placebo'}        
%         subidx=find(sublist.(group{1})==excludesubs(i));
%         if subidx
%             for j=1:length(variablenames)
%                 R_vol_adpt_othervol_model.(variablenames{j}).(group{1})(:,subidx,:)=nan;
%             end
%         end
%         clear subidx
%     end
% end

%stats
[h,p]=ttest2(R_vol_adpt_othervol_model.win_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.win_vol_adpt_tr.placebo(1,:),...
    R_vol_adpt_othervol_model.win_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.win_vol_adpt_tr.reboxetine(1,:))
[h,p]=ttest2(R_vol_adpt_othervol_model.win_sta_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.win_sta_adpt_tr.placebo(1,:),...
    R_vol_adpt_othervol_model.win_sta_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.win_sta_adpt_tr.reboxetine(1,:))

[h,p]=ttest2(R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(1,:),...
    R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(1,:))
[h,p]=ttest2(R_vol_adpt_othervol_model.loss_sta_adpt_tr.placebo(2,:)-R_vol_adpt_othervol_model.loss_sta_adpt_tr.placebo(1,:),...
    R_vol_adpt_othervol_model.loss_sta_adpt_tr.reboxetine(2,:)-R_vol_adpt_othervol_model.loss_sta_adpt_tr.reboxetine(1,:))

[h,p]=ttest2(R_vol_adpt_othervol_model.betas.placebo(2,:,3)-R_vol_adpt_othervol_model.betas.placebo(1,:,3),...
    R_vol_adpt_othervol_model.betas.reboxetine(2,:,3)-R_vol_adpt_othervol_model.betas.reboxetine(1,:,3))

[h,p]=ttest2(log(R_vol_adpt_othervol_model.betas.placebo(2,:,3))-log(R_vol_adpt_othervol_model.betas.placebo(1,:,3)),...
    log(R_vol_adpt_othervol_model.betas.reboxetine(2,:,3))-log(R_vol_adpt_othervol_model.betas.reboxetine(1,:,3)))

[h,p]=ttest2(R_vol_adpt_othervol_model.loss_vol_adpt_tr.placebo(2,:),...
    R_vol_adpt_othervol_model.loss_vol_adpt_tr.reboxetine(2,:))
[h,p]=ttest(R_vol_adpt_othervol_model.loss_sta_adpt_tr.placebo(1,:))
plt_2lr(R_vol_adpt_othervol_model.winalphas,R_vol_adpt_othervol_model.lossalphas,...
    {'placebo'},1,'vol_adpt_by_othervol',figdir);
plt_2lr(R_vol_adpt_othervol_model.winalphas,R_vol_adpt_othervol_model.lossalphas,...
    {'placebo'},2,'vol_adpt_by_othervol',figdir);
plt_2lr(R_vol_adpt_othervol_model.winalphas,R_vol_adpt_othervol_model.lossalphas,...
    {'reboxetine'},1,'vol_adpt_by_othervol',figdir);
plt_2lr(R_vol_adpt_othervol_model.winalphas,R_vol_adpt_othervol_model.lossalphas,...
    {'reboxetine'},2,'vol_adpt_by_othervol',figdir);
plt_r_lrs_adapt_to_othervol(R_vol_adpt_othervol_model,{'reboxetine'},1,'vol_adpt_by_othervol',figdir)
plt_r_lrs_adapt_to_othervol(R_vol_adpt_othervol_model,{'reboxetine'},2,'vol_adpt_by_othervol',figdir)
plt_r_lrs_adapt_to_othervol(R_vol_adpt_othervol_model,{'placebo'},1,'vol_adpt_by_othervol',figdir)
plt_r_lrs_adapt_to_othervol(R_vol_adpt_othervol_model,{'placebo'},2,'vol_adpt_by_othervol',figdir)

plt_r_lrs_adapt_to_othervol_diff(R_vol_adpt_othervol_model,{'placebo'},'vol_adpt_by_othervol',figdir)
plt_r_lrs_adapt_to_othervol_diff(R_vol_adpt_othervol_model,{'reboxetine'},'vol_adpt_by_othervol',figdir)

plt_barplot_lr_adaptation_reb_pla
plt_1beta(R_vol_adpt_othervol_model.betas,{'reboxetine'},1,'vol_adapt_model',figdir)
plt_1beta(R_vol_adpt_othervol_model.betas,{'reboxetine'},2,'vol_adapt_model',figdir)
plt_1beta(R_vol_adpt_othervol_model.betas,{'placebo'},1,'vol_adapt_model',figdir)
plt_1beta(R_vol_adpt_othervol_model.betas,{'placebo'},2,'vol_adapt_model',figdir)
% write results for spss
T_s=array2table(sublist.all,'VariableNames',{'subject'});

T_other_vol_adpt_model_inv=array2table([inv_logit(squeeze(R_vol_adpt_othervol_model.winalphas.all(1,:,:))),...
                                      inv_logit(squeeze(R_vol_adpt_othervol_model.winalphas.all(2,:,:))),...
                                      inv_logit(squeeze(R_vol_adpt_othervol_model.lossalphas.all(1,:,:))),...
                                      inv_logit(squeeze(R_vol_adpt_othervol_model.lossalphas.all(2,:,:))),...
  R_vol_adpt_othervol_model.win_vol_adpt_tr.all(1,:)',R_vol_adpt_othervol_model.win_vol_adpt_tr.all(2,:)',...
  R_vol_adpt_othervol_model.win_sta_adpt_tr.all(1,:)',R_vol_adpt_othervol_model.win_sta_adpt_tr.all(2,:)',...
  R_vol_adpt_othervol_model.loss_vol_adpt_tr.all(1,:)',R_vol_adpt_othervol_model.loss_vol_adpt_tr.all(2,:)',...
  R_vol_adpt_othervol_model.loss_sta_adpt_tr.all(1,:)',R_vol_adpt_othervol_model.loss_sta_adpt_tr.all(2,:)',...
  log(squeeze(R_vol_adpt_othervol_model.betas.all(1,:,:))),log(squeeze(R_vol_adpt_othervol_model.betas.all(2,:,:)))],...
    'VariableNames',{'rew_alpha_inv_bothv_visit1','rew_alpha_inv_winv_visit1','rew_alpha_inv_lossv_visit1','rew_alpha_inv_boths_visit1',...
                                  'rew_alpha_inv_bothv_visit2','rew_alpha_inv_winv_visit2','rew_alpha_inv_lossv_visit2','rew_alpha_inv_boths_visit2',...
                                  'loss_alpha_inv_bothv_visit1','loss_alpha_inv_winv_visit1','loss_alpha_inv_lossv_visit1','loss_alpha_inv_boths_visit1',...
                                  'loss_alpha_inv_bothv_visit2','loss_alpha_inv_winv_visit2','loss_alpha_inv_lossv_visit2','loss_alpha_inv_boths_visit2',...
                                  'rew_alpha_adpt_lossv_visit1','rew_alpha_adpt_lossv_visit2',...
                                  'rew_alpha_adpt_losss_visit1','rew_alpha_adpt_losss_visit2',...
                                  'loss_alpha_adpt_winv_visit1','loss_alpha_adpt_winv_visit2',...
                                  'loss_alpha_adpt_wins_visit1','loss_alpha_adpt_wins_visit2'...
                                  'beta_log_bothv_visit1','beta_log_winv_visit1','beta_log_lossv_visit1','beta_log_boths_visit1',...
                                  'beta_log_bothv_visit2','beta_log_winv_visit2','beta_log_lossv_visit2','beta_log_boths_visit2'});
T_ques=array2table([Ques.treatment,Ques.QIDS_v1,Ques.QIDS_v2,Ques.STAI_S_v1,Ques.STAI_S_v2,...
    Ques.STAI_T_v1,Ques.STAI_T_v2],...
    'VariableNames',{'treatment','QIDS_v1','QIDS_v2','sSTAI_v1','sSTAI_v2','tSTAI_v1','tSTAI_v2'});
T=[T_s,T_ques,T_other_vol_adpt_model_inv];
writetable(T,[datadir,'2019_reboxetine_othervol_adpt_4betas_model_results_alpha_std2_beta_std1'])
%% model fitting for all data in data filefolder
for group={'reboxetine','placebo'}
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
%sort the information and choice made for each block condition
information=[data.winpos(((blkindex(i)-1)*tn+1):blkindex(i)*tn),data.losspos(((blkindex(i)-1)*tn+1):blkindex(i)*tn)];
choice=data.choice(((blkindex(i)-1)*tn+1):blkindex(i)*tn);

%fit models
%2alpha2beta
result_2lr_2b.(group{1})(v,ss,i,:)=fit_linked_2lr_2beta_il_rev(information,choice, start, abandontn);
%1alpha1beta
result_1lr_1b.(group{1})(v,ss,i,:)=fit_linked_1lr_1beta_il(information,choice, start, abandontn);
%2alpha1beta
result_2lr_1b.(group{1})(v,ss,i,:)=fit_linked_2lr_1beta_il(information,choice, start, abandontn);
%2alphas1beta+1bias term
result_2lr_1b_bias.(group{1})(v,ss,i,:)=fit_linked_2lr_beta_add(information,choice, start, abandontn);
%2alpha1beta+choice_stickness
result_2lr_1b_choice_stickness.(group{1})(v,ss,i,:)=fit_linked_2lr_1beta_il_choice_stickness(information,choice, start, abandontn);
%opt1result positive negative PE learning rates
result_PNPE_2lr_1b_opt1.(group{1})(v,ss,i,:)=fit_linked_PNPE_2lr_1beta_opt1results(information,choice, 0, abandontn);
%model with seperate alphas for positive and negative prediction errors
result_PNPE_2lr_2b.(group{1})(v,ss,i,:)=fit_linked_PNPE_2lr_2beta_chosen(information,choice, start, abandontn);
end
end
end
end
%%
save('model_results.mat')
%load('model_results.mat')
%% plot bargraph for mean alpha and beta for 2lr 2beta model
for group={'reboxetine','placebo'}
    for v=1:size(result_2lr_2b.(group{1}),1)
        for i=1:size(result_2lr_2b.(group{1}),2)
            for j=1:size(result_2lr_2b.(group{1}),3)
                  alpha_rew_2lr_2b.(group{1})(v,i,j)=getfield(result_2lr_2b.(group{1}),{v,i,j},'mean_alpha_rew');
                  alpha_loss_2lr_2b.(group{1})(v,i,j)=getfield(result_2lr_2b.(group{1}),{v,i,j},'mean_alpha_loss');
                  beta_rew_2lr_2b.(group{1})(v,i,j)=getfield(result_2lr_2b.(group{1}),{v,i,j},'mean_beta_rew');
                  beta_loss_2lr_2b.(group{1})(v,i,j)=getfield(result_2lr_2b.(group{1}),{v,i,j},'mean_beta_loss'); 
            end
        end
    end
end
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_2lr_2b
        plot_bargarph_betas_rescorla_wagner_2lr_2b
    end
    plot_bargarph_alphas_rescorla_wagner_2lr_2b_dif
    plot_bargarph_betas_rescorla_wagner_2lr_2b_diff
end

%% plot bargraph for mean alpha and beta for 2lr 1beta model
for group={'reboxetine','placebo'}
    for v=1:size(result_2lr_1b.(group{1}),1)
        for i=1:size(result_2lr_1b.(group{1}),2)
            for j=1:size(result_2lr_1b.(group{1}),3)
                  alpha_rew_2lr_1b.(group{1})(v,i,j)=getfield(result_2lr_1b.(group{1}),{v,i,j},'mean_alpha_rew');
                  alpha_loss_2lr_1b.(group{1})(v,i,j)=getfield(result_2lr_1b.(group{1}),{v,i,j},'mean_alpha_loss'); 
                  beta_2lr_1b.(group{1})(v,i,j)=getfield(result_2lr_1b.(group{1}),{v,i,j},'mean_beta'); 
            end
        end
    end
end
%
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_2lr_1b
        plot_bargarph_betas_rescorla_wagner_2lr_1b
        plot_bargarph_diff_alphas_rescorla_wagner_2lr_1b
    end
    plot_bargarph_alphas_rescorla_wagner_2lr_1b_dif
    plot_bargarph_betas_rescorla_wagner_2lr_1b_diff
    plot_bargarph_diff_alphas_rescorla_wagner_2lr_1b_diff
end
%% plot bargraph for mean alpha and beta for 2lr 1beta + bias model
for group={'reboxetine','placebo'}
    for v=1:size(result_2lr_1b_bias.(group{1}),1)
        for i=1:size(result_2lr_1b_bias.(group{1}),2)
            for j=1:size(result_2lr_1b_bias.(group{1}),3)
                  alpha_rew_2lr_1b_bias.(group{1})(v,i,j)=getfield(result_2lr_1b_bias.(group{1}),{v,i,j},'mean_alpha_rew');
                  alpha_loss_2lr_1b_bias.(group{1})(v,i,j)=getfield(result_2lr_1b_bias.(group{1}),{v,i,j},'mean_alpha_loss'); 
                  beta_2lr_1b_bias.(group{1})(v,i,j)=getfield(result_2lr_1b_bias.(group{1}),{v,i,j},'mean_beta'); 
                  bias_2lr_1b_bias.(group{1})(v,i,j)=getfield(result_2lr_1b_bias.(group{1}),{v,i,j},'mean_val_add'); 
            end
        end
    end
end
%
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_2lr_1b_bias
        plot_bargarph_betas_rescorla_wagner_2lr_1b_bias
        plot_bargarph_biasterm_rescorla_wagner_2lr_1b_bias
        plot_bargarph_diff_alphas_rescorla_wagner_2lr_1b_bias
    end
    plot_bargarph_alphas_rescorla_wagner_2lr_1b_bias_dif
   plot_bargarph_betas_rescorla_wagner_2lr_1b_bias_diff
   plot_bargarph_biasterm_rescorla_wagner_2lr_1b_bias_diff
   plot_bargarph_diff_alphas_rescorla_wagner_2lr_1b_bias_diff
end
%% plot bargraph for mean alpha and beta for 2lr 1beta model + choice stickness
for group={'reboxetine','placebo'}
    for v=1:size(result_2lr_1b_choice_stickness.(group{1}),1)
        for i=1:size(result_2lr_1b_choice_stickness.(group{1}),2)
            for j=1:size(result_2lr_1b_choice_stickness.(group{1}),3)
                  alpha_rew_2lr_1b_choice_stickness.(group{1})(v,i,j)=getfield(result_2lr_1b_choice_stickness.(group{1}),{v,i,j},'mean_alpha_rew');
                  alpha_loss_2lr_1b_choice_stickness.(group{1})(v,i,j)=getfield(result_2lr_1b_choice_stickness.(group{1}),{v,i,j},'mean_alpha_loss'); 
                  beta_2lr_1b_choice_stickness.(group{1})(v,i,j)=getfield(result_2lr_1b_choice_stickness.(group{1}),{v,i,j},'mean_beta'); 
                  stk_2lr_1b_choice_stickness.(group{1})(v,i,j)=getfield(result_2lr_1b_choice_stickness.(group{1}),{v,i,j},'mean_stickness'); 
            end
        end
    end
end
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_2lr_1b_choice_stickness
        plot_bargarph_betas_rescorla_wagner_2lr_1b_choice_stickness
        plot_bargarph_stickness_rescorla_wagner_2lr_1b_choice_stickness
    end
     plot_bargarph_alphas_rescorla_wagner_2lr_1b_choice_stk_diff
     plot_bargarph_betas_rescorla_wagner_2lr_1b_choice_stk_diff
     plot_bargarph_stickness_rescorla_wagner_2lr_1b_choice_stk_diff
end

%% 1lr1b
for group={'reboxetine','placebo'}
    for v=1:size(result_1lr_1b.(group{1}),1)
        for i=1:size(result_1lr_1b.(group{1}),2)
            for j=1:size(result_1lr_1b.(group{1}),3)
                  alpha_1lr_1b.(group{1})(v,i,j)=getfield(result_1lr_1b.(group{1}),{v,i,j},'mean_alpha'); 
                  beta_1lr_1b.(group{1})(v,i,j)=getfield(result_1lr_1b.(group{1}),{v,i,j},'mean_beta'); 
            end
        end
    end
end
%
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_1lr_1b
        plot_bargarph_betas_rescorla_wagner_1lr_1b
    end
    plot_bargarph_alphas_rescorla_wagner_1lr_1b_diff
    plot_bargarph_betas_rescorla_wagner_1lr_1b_diff
end
%% PNPE 2lr1b opt1
for group={'reboxetine','placebo'}
    for v=1:size(result_PNPE_2lr_1b_opt1.(group{1}),1)
        for i=1:size(result_PNPE_2lr_1b_opt1.(group{1}),2)
            for j=1:size(result_PNPE_2lr_1b_opt1.(group{1}),3)
                  alpha_PPE_2lr_1b_opt1.(group{1})(v,i,j)=getfield(result_PNPE_2lr_1b_opt1.(group{1}),{v,i,j},'mean_alpha_pos');
                  alpha_NPE_2lr_1b_opt1.(group{1})(v,i,j)=getfield(result_PNPE_2lr_1b_opt1.(group{1}),{v,i,j},'mean_alpha_neg');
                  beta_PNPE_2lr_1b_opt1.(group{1})(v,i,j)=getfield(result_PNPE_2lr_1b_opt1.(group{1}),{v,i,j},'mean_beta');
            end
        end
    end
end
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_opt1_PNPE_2lr_1b
        plot_bargarph_betas_rescorla_wagner_opt1_PNPE_2lr_1b
    end
    plot_bargarph_alphas_rescorla_wagner_opt1_PNPE_2lr_1b_diff
    plot_bargarph_betas_rescorla_wagner_opt1_PNPE_2lr_1b_diff
end
%% separete alphas for negative positive Prediction errors
for group={'reboxetine','placebo'}
    for v=1:size(result_PNPE_2lr_2b.(group{1}),1)
        for i=1:size(result_PNPE_2lr_2b.(group{1}),2)
            for j=1:size(result_PNPE_2lr_2b.(group{1}),3)
                  alpha_PPE_2lr_2b.(group{1})(v,i,j)=getfield(result_PNPE_2lr_2b.(group{1}),{v,i,j},'mean_alpha_pos');
                  alpha_NPE_2lr_2b.(group{1})(v,i,j)=getfield(result_PNPE_2lr_2b.(group{1}),{v,i,j},'mean_alpha_neg');
                  beta_rew_PNPE_2lr_2b.(group{1})(v,i,j)=getfield(result_PNPE_2lr_2b.(group{1}),{v,i,j},'mean_beta_rew');
                  beta_loss_PNPE_2lr_2b.(group{1})(v,i,j)=getfield(result_PNPE_2lr_2b.(group{1}),{v,i,j},'mean_beta_loss'); 
            end
        end
    end
end
for group={'reboxetine','placebo'}
    for v=1:2
        plot_bargarph_alphas_rescorla_wagner_PNPE_2lr_2b
        plot_bargarph_betas_rescorla_wagner_PNPE_2lr_2b
    end
    plot_bargarph_alphas_rescorla_wagner_PNPE_2lr_2b_diff
    plot_bargarph_betas_rescorla_wagner_PNPE_2lr_2b_diff
end
%% compare mean BICs for all models
for group={'reboxetine','placebo'}
    for v=1:size(result_2lr_2b.(group{1}),1)
        for i=1:size(result_2lr_2b.(group{1}),2)
            for j=1:size(result_2lr_2b.(group{1}),3)
                  BIC_2lr_2b.(group{1})(v,i,j)=getfield(result_2lr_2b.(group{1}),{v,i,j},'BIC');
                  BIC_2lr_1b.(group{1})(v,i,j)=getfield(result_2lr_1b.(group{1}),{v,i,j},'BIC');
                  BIC_1lr_1b.(group{1})(v,i,j)=getfield(result_1lr_1b.(group{1}),{v,i,j},'BIC');
                  BIC_PNPE_2lr_2b.(group{1})(v,i,j)=getfield(result_PNPE_2lr_2b.(group{1}),{v,i,j},'BIC');
                  BIC_PNPE_2lr_1b_opt1.(group{1})(v,i,j)=getfield(result_PNPE_2lr_1b_opt1.(group{1}),{v,i,j},'BIC');
                  BIC_2lr_1b_choice_stickness.(group{1})(v,i,j)=getfield(result_2lr_1b_choice_stickness.(group{1}),{v,i,j},'BIC');
                  BIC_2lr_1b_bias.(group{1})(v,i,j)=getfield(result_2lr_1b_bias.(group{1}),{v,i,j},'BIC');
            end
        end
    end
end
BIC_2lr_2b.all=[BIC_2lr_2b.reboxetine,BIC_2lr_2b.placebo];
BIC_2lr_1b.all=[BIC_2lr_1b.reboxetine,BIC_2lr_1b.placebo];
BIC_1lr_1b.all=[BIC_1lr_1b.reboxetine,BIC_1lr_1b.placebo];
BIC_PNPE_2lr_2b.all=[BIC_PNPE_2lr_2b.reboxetine,BIC_PNPE_2lr_2b.placebo];
BIC_PNPE_2lr_1b_opt1.all=[BIC_PNPE_2lr_1b_opt1.reboxetine,BIC_PNPE_2lr_1b_opt1.placebo];
BIC_2lr_1b_choice_stickness.all=[BIC_2lr_1b_choice_stickness.reboxetine,BIC_2lr_1b_choice_stickness.placebo];
BIC_2lr_1b_bias.all=[BIC_2lr_1b_bias.reboxetine,BIC_2lr_1b_bias.placebo];
for group={'reboxetine','placebo','all'}
    for v=1:2
plot_bargarph_BICs
    end
end
%% AIC
for i=1:size(result_2lr_2b,1)
    for j=1:4%size(result_2lr_2b,2)
          AIC_2lr_2b(v,i,j)=getfield(result_2lr_2b,{i,j+1},'AIC');
          AIC_2lr_1b(v,i,j)=getfield(result_2lr_1b,{i,j+1},'AIC');
          AIC_1lr_1b(v,i,j)=getfield(result_1lr_1b,{i,j+1},'AIC');
          AIC_1lr_1b_opt1(v,i,j)=getfield(result_PNPE_2lr_1b_opt1,{i,j+1},'AIC');
          AIC_PNPE_2lr_2b(v,i,j)=getfield(result_PNPE_2lr_2b,{i,j+1},'AIC');
          AIC_PNPE_2lr_1b(v,i,j)=getfield(result_PNPE_2lr_1b,{i,j+1},'AIC');
          AIC_PNPE_2lr_1b_opt1(v,i,j)=getfield(result_PNPE_2lr_1b_opt1,{i,j+1},'AIC');
          AIC_2lr_1b_choice_stickness(v,i,j)=getfield(result_2lr_1b_choice_stickness,{i,j+1},'AIC');
    end
end
plot_bargarph_AICs

%%
T_s=table(Ques.subnum,Ques.treatment,'VariableNames',{'subjectnum','treatment'});
T_rew_alphas=array2table([squeeze(inv_logit(alpha_rew_2lr_1b.all(1,:,:))),squeeze(inv_logit(alpha_rew_2lr_1b.all(2,:,:)))],...
    'VariableNames',{'rew_alpha_2lr1b_visit1_bothv','rew_alpha_2lr1b_visit1_winv','rew_alpha_2lr1b_visit1_lossv','rew_alpha_2lr1b_visit1_boths',...
                                    'rew_alpha_2lr1b_visit2_bothv','rew_alpha_2lr1b_visit2_winv','rew_alpha_2lr1b_visit2_lossv','rew_alpha_2lr1b_visit2_boths'});
T_loss_alphas=array2table([squeeze(inv_logit(alpha_loss_2lr_1b.all(1,:,:))),squeeze(inv_logit(alpha_loss_2lr_1b.all(2,:,:)))],...
    'VariableNames',{'loss_alpha_2lr1b_visit1_bothv','loss_alpha_2lr1b_visit1_winv','loss_alpha_2lr1b_visit1_lossv','loss_alpha_2lr1b_visit1_boths',...
                                    'loss_alpha_2lr1b_visit2_bothv','loss_alpha_2lr1b_visit2_winv','loss_alpha_2lr1b_visit2_lossv','loss_alpha_2lr1b_visit2_boths'});
T_betas=array2table([squeeze(log(beta_2lr_1b.all(1,:,:))),squeeze(log(beta_2lr_1b.all(2,:,:)))],...
    'VariableNames',{'beta_2lr1b_visit1_bothv','beta_2lr1b_visit1_winv','beta_2lr1b_visit1_lossv','beta_2lr1b_visit1_boths',...
                                    'beta_2lr1b_visit2_bothv','beta_2lr1b_visit2_winv','beta_2lr1b_visit2_lossv','beta_2lr1b_visit2_boths'});  
%rm anova model for alphas
T_stat_alphas=[T_s,T_rew_alphas,T_loss_alphas];
within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
     'VariableNames',{'valence','visit','winvol','lossvol'});
 rm=fitrm(T_stat_alphas,'rew_alpha_2lr1b_visit1_bothv-loss_alpha_2lr1b_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl_alpha_2lr_1b=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
%rm anova model for batas
T_stat_betas=[T_s,T_betas];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
     'VariableNames',{'visit','winvol','lossvol'});
rm=fitrm(T_stat_betas,'beta_2lr1b_visit1_bothv-beta_2lr1b_visit2_boths ~ treatment','WithinDesign',within);
ranovatbl_beta_2lr1b=ranova(rm,'WithinModel','visit*winvol*lossvol')
T=[T_s, T_rew_alphas,T_loss_alphas,T_betas];
 writetable(T,[datadir,'2lr1b_model'])
%% do stats for 2lr1b bias model
T_s=table(Ques.subnum,Ques.treatment,'VariableNames',{'subjectnum','treatment'});
T_rew_alphas=array2table([squeeze(inv_logit(alpha_rew_2lr_1b_bias.all(1,:,:))),squeeze(inv_logit(alpha_rew_2lr_1b_bias.all(2,:,:)))],...
    'VariableNames',{'rew_alpha_2lr1b_bias_visit1_bothv','rew_alpha_2lr1b_bias_visit1_winv','rew_alpha_2lr1b_bias_visit1_lossv','rew_alpha_2lr1b_bias_visit1_boths',...
                                    'rew_alpha_2lr1b_bias_visit2_bothv','rew_alpha_2lr1b_bias_visit2_winv','rew_alpha_2lr1b_bias_visit2_lossv','rew_alpha_2lr1b_bias_visit2_boths'});
T_loss_alphas=array2table([squeeze(inv_logit(alpha_loss_2lr_1b_bias.all(1,:,:))),squeeze(inv_logit(alpha_loss_2lr_1b_bias.all(2,:,:)))],...
    'VariableNames',{'loss_alpha_2lr1b_bias_visit1_bothv','loss_alpha_2lr1b_bias_visit1_winv','loss_alpha_2lr1b_bias_visit1_lossv','loss_alpha_2lr1b_bias_visit1_boths',...
                                    'loss_alpha_2lr1b_bias_visit2_bothv','loss_alpha_2lr1b_bias_visit2_winv','loss_alpha_2lr1b_bias_visit2_lossv','loss_alpha_2lr1b_bias_visit2_boths'});
T_betas=array2table([squeeze(log(beta_2lr_1b_bias.all(1,:,:))),squeeze(log(beta_2lr_1b_bias.all(2,:,:)))],...
    'VariableNames',{'beta_2lr1b_bias_visit1_bothv','beta_2lr1b_bias_visit1_winv','beta_2lr1b_bias_visit1_lossv','beta_2lr1b_bias_visit1_boths',...
                                    'beta_2lr1b_bias_visit2_bothv','beta_2lr1b_bias_visit2_winv','beta_2lr1b_bias_visit2_lossv','beta_2lr1b_bias_visit2_boths'});
T_biases=array2table([squeeze(abs(bias_2lr_1b_bias.all(1,:,:))),squeeze(abs(bias_2lr_1b_bias.all(2,:,:)))],...
    'VariableNames',{'bias_2lr1b_bias_visit1_bothv','bias_2lr1b_bias_visit1_winv','bias_2lr1b_bias_visit1_lossv','bias_2lr1b_bias_visit1_boths',...
                                    'bias_2lr1b_bias_visit2_bothv','bias_2lr1b_bias_visit2_winv','bias_2lr1b_bias_visit2_lossv','bias_2lr1b_bias_visit2_boths'}); 
                                
%rm anova model for alphas
T_stat_alphas=[T_s,T_rew_alphas,T_loss_alphas];
within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
     'VariableNames',{'valence','visit','winvol','lossvol'});
 rm=fitrm(T_stat_alphas,'rew_alpha_2lr1b_bias_visit1_bothv-loss_alpha_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl_alphas_2lr1b_bias=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
%rm anova model for batas
T_stat_betas=[T_s,T_betas];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
     'VariableNames',{'visit','winvol','lossvol'});
rm=fitrm(T_stat_betas,'beta_2lr1b_bias_visit1_bothv-beta_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
ranovatbl_betas_2lr1b_bias=ranova(rm,'WithinModel','visit*winvol*lossvol')
%rm anova model for bias term
T_stat_biases=[T_s,T_biases];
within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
     'VariableNames',{'visit','winvol','lossvol'});
 rm=fitrm(T_stat_biases,'bias_2lr1b_bias_visit1_bothv-bias_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
 ranovatbl_bias_2lr1b_bias=ranova(rm,'WithinModel','visit*winvol*lossvol')
 
 % write table to txt for spss
 T=[T_s, T_rew_alphas,T_loss_alphas,T_betas,T_biases];
 writetable(T,[datadir,'2lr1b_bias_model'])