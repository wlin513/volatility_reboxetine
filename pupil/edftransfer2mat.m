clear 
clc
getfolders;
%%
addpath([codedir,'/edf-converter-master'])
%%
filelist=dir([datadir,'*/*.edf']);

for s=1:length(filelist)
    dat=Edf2Mat([filelist(s).folder,'/',filelist(s).name]);
    display(filelist(s).name);
    save(dat);
    clear dat;
end
disp('complete');