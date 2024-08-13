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

        sub_data.tracker_data=out;
        store=parsemisstracker(sub_data);


        for i=1:size(store.options_missa,1)    
            rewards_missing(sub,i)=sum(store.rew_outcome_missb(i,:))/size(store.options_missa,2)*100;
            loss_missing(sub,i)=sum(store.pun_outcome_missb(i,:))/size(store.options_missa,2)*100;   
        end

        rew_exclude(1,320)=zeros;
        loss_exclude(1,320)=zeros;
        for k=1:size(store.options_missa,1) 

            if rewards_missing(sub,k)<=50
                rew_exclude(k)=0;
            else
                rew_exclude(k)=1;
            end

            if loss_missing(sub,k)<=50
                pun_exclude(k)=0;
            else
                pun_exclude(k)=1;
            end
        end

        save([datadir,name,'_visit_',num2str(visit), '_pupil_exclude_trials.mat'],'rew_exclude','pun_exclude')
    end
end

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
% sub_data.tracker_data=out;
% store=parsemisstracker(sub_data);
% 
% 
% for i=1:size(store.options_missa,1)    
%     rewards_missing(sub,i)=sum(store.rew_outcome_missb(i,:))/size(store.options_missa,2)*100;
%     loss_missing(sub,i)=sum(store.pun_outcome_missb(i,:))/size(store.options_missa,2)*100;   
% end
% 
% 
% rew_exclude=[];
% loss_exclude=[];
% for k=1:size(store.options_missa,1) 
%     
%     if rewards_missing(sub,k)<=50
%         rew_exclude(k)=0;
%     else
%         rew_exclude(k)=1;
%     end
%     
%     if loss_missing(sub,k)<=50
%         pun_exclude(k)=0;
%     else
%         pun_exclude(k)=1;
%     end
% end
% name=subs{sub};
% dir_name=['C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\',name,'\'];
% 
% pth=dir_name;
% cd(pth)
% 
% save([name '_pupil_exclude_trials.mat'],'rew_exclude','pun_exclude')
% 
% end
% 
