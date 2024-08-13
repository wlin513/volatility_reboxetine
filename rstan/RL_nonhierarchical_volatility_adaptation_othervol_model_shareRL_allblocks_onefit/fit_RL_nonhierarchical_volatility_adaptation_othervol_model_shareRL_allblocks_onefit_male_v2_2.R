# script to read in data from 2018_fMRI single subject data and fit parameters

rm(list=ls())


# configure filepaths -----------------------------------------------------

setwd("/home/wlin/Documents/2019_drug_study/code/rstan")
#where the data files are
data_dir<-"/home/wlin/Documents/2019_drug_study/data"
dest_dir<-"/home/wlin/Documents/2019_drug_study/data/rstan/RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit"
stan_dir<-"/home/wlin/Documents/2019_drug_study/code/rstan/RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit"
fig_dir<-"/home/wlin/Documents/2019_drug_study/tmp_fig/rstan/RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit"

ifelse(!dir.exists(dest_dir), dir.create(dest_dir), NA)
ifelse(!dir.exists(fig_dir), dir.create(fig_dir),NA)
#some user functions
source(file.path(stan_dir,"fun_read_original_data_files.R"))
source(file.path(stan_dir,"fun_parse_data_into_blocks.R"))
source(file.path(stan_dir,"inv_logit.R"))
source(file.path(stan_dir,"logit.R"))
source(file.path(stan_dir,"read_ests_vol_adpt_othervol_model_shareRL_allblocks.R"))

#libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# configure task details --------------------------------------------------

#num trials in block
ntrials=80
nabandon=10
nblk=4
#the block order information
blockorders<- matrix(c(1,2,3,4, 1,3,2,4,4,2,3,1,4,3,2,1,2,1,4,3,2,4,1,3,3,4,1,2,3,1,4,2),nrow=4,ncol=8)
#read in the list of subjects
subjects<-read.table(file=file.path(data_dir,"rstan_sublist_male_v2_2.txt"))

#total number of subjects
nsubs=nrow(subjects)


#create arrays for data. There are 80 trials per block and 4 blocks

win_outcome<-array(dim=c(ntrials,4,nsubs))
choice<-array(dim=c(ntrials,4,nsubs))
loss_outcome<-array(dim=c(ntrials,4,nsubs))
includeTrial<-array(dim=c(ntrials,4,nsubs))
session<-array(dim=c(nsubs))

# loop through subjects getting their data ----------------------------------------

for(sn in 1:nsubs)
{
  sub_data<-fun_parse_data_into_blocks(fun_read_original_data_files(file.path(data_dir,subjects[sn,1]),sep=""))
  
  #sort data 
  
  blockorder=blockorders[,sub_data$blktype]
  
  i=1
  for(block in blockorder)
  { print(block)
    print(i)
    win_outcome[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_wins",sep="")])
    choice[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_choice",sep="")])
    loss_outcome[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_loss",sep="")])
    ifmakechoice=is.finite(data.matrix(sub_data$data[paste("block",i,"_RT",sep="")]))
    includet=as.numeric(ifmakechoice)
    includet[1:nabandon]=0
    includeTrial[,block,sn]<-includet
    i=i+1
  }
  session[sn]<-sub_data$nvisit
}
# fit stan ----------------------------------------------------------------



#scan blocks
stan_data<-list(ntr=ntrials,nblk=nblk, nsub=nsubs, opt1Chosen=choice, opt1winout=win_outcome, opt1lossout=loss_outcome,includeTrial=includeTrial)
fit_scandata <- stan(file = file.path(stan_dir,'RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit.stan'), data = stan_data, iter = 5000, chains = 4, control=list(adapt_delta=0.99, max_treedepth = 15))
save(fit_scandata,file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit_male_v2_2.RData"))
#library("shinystan")
#my_sso <- launch_shinystan(fit_scandata)



# figures -----------------------------------------------------------------
#read data and save as mat file
#load(file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit.RData"))

IDs=sub('\\/.*',"",subjects$V1)
blocknames=c("both volatile","win volatile","loss volatile","both stable")
alldata=read_ests(fit_scandata,IDs,blocknames)

library(ggplot2)
library(plyr)
library(reshape2)
betas=subset(alldata,alldata$variables=="beta")
win_vol_adpt_tr=subset(alldata,alldata$variables=="win_vol_adpt_tr")
loss_vol_adpt_tr=subset(alldata,alldata$variables=="loss_vol_adpt_tr")
win_sta_adpt_tr=subset(alldata,alldata$variables=="win_sta_adpt_tr")
loss_sta_adpt_tr=subset(alldata,alldata$variables=="loss_sta_adpt_tr")
winalphas=subset(alldata,alldata$variables=="winalpha")
lossalphas=subset(alldata,alldata$variables=="lossalpha")
allalphas=rbind(winalphas,lossalphas)

lgi_alphas=allalphas
lgi_alphas$values=logit(allalphas$values)


alphameans <- ddply(lgi_alphas, c("block", "variables"), summarise,
               mean=inv_logit(mean(values)),se=inv_logit(mean(values)+sd(values)/sqrt(length(values)))-inv_logit(mean(values)))
palpha<-ggplot(alphameans,aes(x=block, y=mean, fill=variables)) +
            xlab("blocks")+ylab("alphas")+
            geom_bar(position=position_dodge(), stat="identity",
            size=.3) +      # Thinner lines
            scale_fill_manual(values=c("lightgreen","red"))+
            theme_minimal()+
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) 

palpha

ggsave("alphas_RL_nonhierarchical_vol_adpt_othervol_model_shareRL_onefit_male_v2.eps", plot = palpha, device = "eps", path = fig_dir,
       scale = 1, width = 20, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)




#library("shinystan")
#my_sso <- launch_shinystan(fit_prescan)
                          
# save data in a mat file -------------------------------------------------
library(R.matlab)
writeMat(file.path(dest_dir,'Result_vol_adpt_othervol_model_shareRL_male_v2_2.mat'),winalphas=winalphas,lossalphas=lossalphas,betas=betas,win_vol_adpt_tr=win_vol_adpt_tr,loss_vol_adpt_tr=loss_vol_adpt_tr,win_sta_adpt_tr=win_sta_adpt_tr,loss_sta_adpt_tr=loss_sta_adpt_tr)
