clear all
getfolders
%%
subs=[1001:1030,2001:2009,2011:2020,2030];
dir_name=[datadir,'ascfiles/'];
for sub=1:size(subs,2)
    tracker_chop_normalise_bet_visits(num2str(subs(sub)),dir_name);
end