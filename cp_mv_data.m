
cd('/home/wlin/Documents/2019_drug_study/code');
clear all
getfolders
%% load subjects and group information
load([datadir,'questionnaire.mat'])
sublist.all=Ques.subnum;
sublist.reboxetine=sublist.all(Ques.treatment==1);
sublist.placebo=sublist.all(Ques.treatment==2);
excludesub=[];
for group={'reboxetine','placebo'}
    for i=1:length(excludesub)
        if isempty(find(sublist.(group{1})==excludesub(i)))
            continue
        else
              excnum(i)=find(sublist.(group{1})==excludesub(i));
        end
    sublist.(group{1})(excnum(i))=[];
    end
end
 group={'placebo'};
 ssublist=sublist.(group{1});
 des_dir=[datadir,'control_behavior_txt_files'];
 mkdir(des_dir)
for sub=1:size(ssublist,1)    
    name=num2str(ssublist(sub));
    sdatadir=[datadir,name,'/'];
    for visit=1:2
         filename=dir([sdatadir,'*_visit_',num2str(visit),'_blktype_*.txt']);
         file=[sdatadir,filename.name];
         copyfile(file,des_dir)
    end
end