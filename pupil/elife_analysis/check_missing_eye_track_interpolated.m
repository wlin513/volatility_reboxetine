%%check the eye tracking data for missing parts
clear all
getfolders
datadir=[datadir,'ascfiles/'];
subs=[1001:1030,2001:2009,2011:2020,2030];
for sub=1:size(subs,2)
         for visit=1:2
            name=num2str(subs(sub));
            file_name=['out_',name,'_visit_',num2str(visit),'.mat'];
            load([datadir,file_name])

        removed=out.removed_data;

        counter=0;
        for i=1:size(removed,1)    
            if removed(i,1)==1 && removed(i,2)==1 
                counter=counter+1;
            end
        end

        missing_percent(sub,visit)=(counter/size(removed,1))*100;
         end
end
% 
% %%check the eye tracking data for missing parts- mag task
% clear all
% subs={'v9','v10','v11','v12','v13','v15','v16','v18', 'v21','v23','v25','v26', 'v28','v29','v31','v34','v38','v40','v42','v43','v44','v45','v47','v48','v51','v52','v54','v55','v56'};
% for sub=1:size(subs,2)
%     
% name=subs{sub};
% dir_name=['C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\',name,'\'];
% 
% pth=dir_name;
% cd(pth)
% load out.mat
% 
% pth='C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\';
% cd(pth)
% 
% removed=out.removed_data;
% 
% counter=0;
% for i=1:size(removed,1)    
%     if removed(i,1)==1 && removed(i,2)==1 
%         counter=counter+1;
%     end
% end
% 
% missing_percent(sub,1)=(counter/size(removed,1))*100;
% end
% 
