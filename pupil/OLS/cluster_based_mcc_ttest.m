function [orisigpoints, bootsigpoints, clustersigpoints] = cluster_based_mcc_ttest(data,sigalpha,b)

%one sample t-test cluster based multiple comparison correction
%
%input should be (1)data(time,sub), (2)significant alpha level, (3)bootstrapping number;
%output would be (1)original significant points(uncorrected student t test), (2)bootstrapping significant points(uncorrected bootstrapping test), 
%(3)significant pointed after cluster based multiple comparision correction


[H,orip,CI,stats] = ttest(data,0,'Dim',2);
orit = stats.tstat;

orisigpoints=(orip<=sigalpha);
orisigpoints=double(orisigpoints);
orisigpoints(orisigpoints==0)=nan;
%% do bootstrap
Nf = size(data,1);

th = b*(1-sigalpha);
% lo = round(b*sigalpha./2);
% hi = b-lo;
% lo = lo+1;

boott = zeros(Nf,b);
bootp = boott;
cdata = data-repmat(mean(data,2),[1 size(data,2)]);

for kk=1:b
    [H,bootp(:,kk),CI,stats] = ttest(cdata(:,randi(size(data,2),size(data,2),1)),0,'Dim',2);
    boott(:,kk) = stats.tstat;
end
%% get bootstrap univariate thresholds

% IMPORTANT: as you will see, on it's own the bootstrap DOES NOT
%            correct for multiple comparisons

sortboot = sort(boott.^2,2); % we use t.^2 (=F) to make our life easier
unibootth = sortboot(:,th); 

bootsigpoints=(orit.^2>=unibootth);
bootsigpoints=double(bootsigpoints);
bootsigpoints(bootsigpoints==0)=nan;
%% get cluster based significant points
% get cluster threshold

clusterbootth = limo_ecluster_make(boott.^2,bootp,sigalpha);

% cluster original data

sigcluster = limo_ecluster_test(orit'.^2,orip',clusterbootth,sigalpha);

% get cluster based significant points
clustersigpoints=sigcluster.elec_mask;
clustersigpoints(clustersigpoints==0)=nan;


