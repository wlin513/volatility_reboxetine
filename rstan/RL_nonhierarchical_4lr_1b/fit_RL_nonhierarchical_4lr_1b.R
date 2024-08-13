# script to read in data from 2018_fMRI single subject data and fit parameters

rm(list=ls())


# configure filepaths -----------------------------------------------------

setwd("/home/wlin/Documents/2019_drug_study/code/rstan")
#where the data files are
data_dir<-"/home/wlin/Documents/2019_drug_study/data"
dest_dir<-"/home/wlin/Documents/2019_drug_study/data/rstan/RL_nonhierarchical_4lr_1b"
stan_dir<-"/home/wlin/Documents/2019_drug_study/code/rstan/RL_nonhierarchical_4lr_1b"
fig_dir<-"/home/wlin/Documents/2019_drug_study/tmp_fig/rstan/RL_nonhierarchical_4lr_1b"

ifelse(!dir.exists(dest_dir), dir.create(dest_dir), NA)
ifelse(!dir.exists(fig_dir), dir.create(fig_dir),NA)
#some user functions
source(file.path(stan_dir,"fun_read_original_data_files.R"))
source(file.path(stan_dir,"fun_parse_data_into_blocks.R"))
source(file.path(stan_dir,"inv_logit.R"))
source(file.path(stan_dir,"logit.R"))
source(file.path(stan_dir,"read_ests_4lr_1b.R"))

#libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# configure task details --------------------------------------------------

#num trials in block
ntrials=80
nabandon=10

#the block order information
blockorders<- matrix(c(1,2,3,4, 1,3,2,4,4,2,3,1,4,3,2,1,2,1,4,3,2,4,1,3,3,4,1,2,3,1,4,2),nrow=4,ncol=8)
#read in the list of subjects
subjects<-read.table(file=file.path(data_dir,"rstan_sublist.txt"))

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

#both vol
stan_data<-list(ntr=ntrials,nsub=nsubs, opt1Chosen=choice[,1,], opt1winout=win_outcome[,1,], opt1lossout=loss_outcome[,1,],includeTrial=includeTrial[,1,])
fit_bothv <- stan(file = file.path(stan_dir,'RL_nonhierarchical_4lr_1b.stan'), data = stan_data, iter = 4000, chains = 4, control=list(adapt_delta=0.999, max_treedepth = 10))
save(fit_bothv,file=file.path(dest_dir,"RL_nonhierarchical_4lr_1b_bothv.RData"))
#library("shinystan")
#my_sso <- launch_shinystan(fit_bothv)


#win vol
stan_data<-list(ntr=ntrials,nsub=nsubs, opt1Chosen=choice[,2,], opt1winout=win_outcome[,2,], opt1lossout=loss_outcome[,2,],includeTrial=includeTrial[,2,])
fit_winv <- stan(file = file.path(stan_dir,'RL_nonhierarchical_4lr_1b.stan'), data = stan_data, iter = 4000, chains = 4, control=list(adapt_delta=0.999, max_treedepth = 10))
save(fit_winv,file=file.path(dest_dir,"RL_nonhierarchical_4lr_1b_winv.RData"))

#library("shinystan")
#my_sso <- launch_shinystan(fit_winv)

#loss vol
stan_data<-list(ntr=ntrials,nsub=nsubs, opt1Chosen=choice[,3,], opt1winout=win_outcome[,3,], opt1lossout=loss_outcome[,3,],includeTrial=includeTrial[,3,])
fit_lossv <- stan(file = file.path(stan_dir,'RL_nonhierarchical_4lr_1b.stan'), data = stan_data, iter = 4000, chains = 4, control=list(adapt_delta=0.999, max_treedepth = 10))
save(fit_lossv,file=file.path(dest_dir,"RL_nonhierarchical_4lr_1b_lossv.RData"))

#library("shinystan")
#my_sso <- launch_shinystan(fit_lossv)

#both stable
stan_data<-list(ntr=ntrials,nsub=nsubs, opt1Chosen=choice[,4,], opt1winout=win_outcome[,4,], opt1lossout=loss_outcome[,4,],includeTrial=includeTrial[,4,])
fit_boths <- stan(file = file.path(stan_dir,'RL_nonhierarchical_4lr_1b.stan'), data = stan_data, iter = 4000, chains = 4, control=list(adapt_delta=0.999, max_treedepth = 10))
save(fit_boths,file=file.path(dest_dir,"RL_nonhierarchical_4lr_1b_boths.RData"))

library("shinystan")
my_sso <- launch_shinystan(fit_boths)

# figures -----------------------------------------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
#load(file=file.path(dest_dir,"RL_nonhierarchical_4lr_1b_lossv.RData"))
IDs=sub('\\/.*',"",subjects$V1)
visit=sub('.*/.*_visit_',"",subjects$V1)
visit=sub('_blktype.*',"",visit)
result_bothv=read_ests(fit_bothv,IDs,"both volatile")
result_winv=read_ests(fit_winv,IDs,"win volatile")
result_lossv=read_ests(fit_lossv,IDs,"loss volatile")
result_boths=read_ests(fit_boths,IDs,"both stable")



alldata<-rbind(result_bothv,result_winv,result_lossv,result_boths)
winposalphas=subset(alldata,alldata$variables=="winposalpha")
winnegalphas=subset(alldata,alldata$variables=="winnegalpha")
lossposalphas=subset(alldata,alldata$variables=="lossposalpha")
lossnegalphas=subset(alldata,alldata$variables=="lossnegalpha")
allalphas=rbind(winposalphas,winnegalphas,lossposalphas,lossnegalphas)
allbetas=subset(alldata,alldata$variables=="beta")
#lgi_alphas=allalphas
#lgi_alphas$values=logit(allalphas$values)
#log_betas=allbetas
#log_betas$values=log(allbetas$values)

#alphameans <- ddply(lgi_alphas, c("block", "variables"), summarise,
#               mean=inv_logit(mean(values)),se=inv_logit(mean(values)+sd(values)/sqrt(length(values)))-inv_logit(mean(values)))
#palpha<-ggplot(alphameans,aes(x=block, y=mean, fill=variables)) +
#            xlab("blocks")+ylab("alphas")+
#            geom_bar(position=position_dodge(), stat="identity",
#            size=.3) +      # Thinner lines
#            scale_fill_manual(values=c("lightgreen","green","lightpink","red"))+
#           theme_minimal()+
#           geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
#                size=.3,    # Thinner lines
#                width=.2,
#               position=position_dodge(.9)) 

#palpha

#ggsave("alphas_RL_nonhierarchical_4lr_1b.eps", plot = palpha, device = "eps", path = fig_dir,
#       scale = 1, width = 20, height = 10, units = "cm",
#      dpi = 300, limitsize = TRUE)

#betameans <- ddply(log_betas, c("block", "variables"), summarise,
#                    mean=exp(mean(values)),se=exp(mean(values)+sd(values)/sqrt(length(values)))-exp(mean(values)))
#pbeta<-ggplot(betameans,aes(x=block, y=mean, fill=variables)) +
#  xlab("blocks")+ylab("betas")+
#  geom_bar(position=position_dodge(), stat="identity",
#           size=.3) +      # Thinner lines
# theme_minimal()+
# geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
#                size=.3,    # Thinner lines
#                width=.2,
#                position=position_dodge(.9)) 

#pbeta

#ggsave("betas_RL_nonhierarchical_4lr_1b.eps", plot = pbeta, device = "eps", path = fig_dir,
#      scale = 1, width = 20, height = 10, units = "cm",
#      dpi = 300, limitsize = TRUE)

#library("shinystan")
#my_sso <- launch_shinystan(fit_prescan)

# save data in a mat file -------------------------------------------------
library(R.matlab)
#load(file.path(dest_dir,"RL_nonhierarchical_4lr_1b_prescan.RData"))
#load(file.path(dest_dir,"RL_nonhierarchical_4lr_1b_bothv.RData"))
#load(file.path(dest_dir,"RL_nonhierarchical_4lr_1b_winv.RData"))
#load(file.path(dest_dir,"RL_nonhierarchical_4lr_1b_lossv.RData"))
#load(file.path(dest_dir,"RL_nonhierarchical_4lr_1b_boths.RData"))
writeMat(file.path(dest_dir,'R_result_4lr_1b.mat'),winposalphas=winposalphas,winnegalphas=winnegalphas,lossposalphas=lossposalphas,lossnegalphas=lossnegalphas)