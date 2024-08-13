fun_read_original_data_files <- function(file.name, ...){
  lines <- scan(file.name, what="character", sep="\n")
  first.line <- min(grep("Trialnumber", lines))
  sub_num<-substr(lines[1],13,16)
  visitnum<-substr(lines[1],28,28)
  blktype<-substr(lines[1],40,40)
  ifreverse<-substr(lines[1],52,52)
  nsubnum=as.numeric(sub_num)
  nvisit=as.numeric(visitnum)
  nblktype=as.numeric(blktype)
  nver=as.numeric(ifreverse)
  ret_list <- list("data"=read.delim(textConnection(lines), skip=first.line-1),"subnum"=nsubnum,"nvisit"=nvisit,"blktype"=nblktype,'ver'=nver)
  return(ret_list)
  }
