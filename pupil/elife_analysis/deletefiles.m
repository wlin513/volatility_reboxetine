getfolders
for i=1:32
target=[datadir, 'elife_analysis/',num2str(i),'/out_tfMRI_vol_',num2str(i),'.mat'];

delete(target);
end
