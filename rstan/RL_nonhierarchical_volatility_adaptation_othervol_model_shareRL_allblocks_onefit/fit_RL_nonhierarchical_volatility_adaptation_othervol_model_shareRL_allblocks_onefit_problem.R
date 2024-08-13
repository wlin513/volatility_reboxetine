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
nabandon=5
nblk=4
#the block order information
blockorders<- matrix(c(1,2,3,4, 1,3,2,4,4,2,3,1,4,3,2,1,2,1,4,3,2,4,1,3,3,4,1,2,3,1,4,2),nrow=4,ncol=8)
#read in the list of subjects
subjects<-read.table(file=file.path(data_dir,"rstan_sublist_problems2.txt"))

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
fit_scandata <- stan(file = file.path(stan_dir,'RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit_problem.stan'), data = stan_data, iter = 20000, chains = 4, control=list(adapt_delta=0.99, max_treedepth = 10))
save(fit_scandata,file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit_problem.RData"))
# dest_dir<-"/home/wlin/Documents/2019_drug_study/data/rstan/RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit"
# load(file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_othervol_model_shareRL_allblocks_onefit_problem.RData"))
# library("shinystan")
# my_sso <- launch_shinystan(fit_scandata)
