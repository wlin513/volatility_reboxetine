function out=OptimalBayesLearner_flat_func(info,flattenpoints)

% version of the Behrens optimal learner which accounts for change in
% stimuli at certain trials. At these points (called flattenpoints) the dimension of
% the learner which codes belief about the stimuli (i.e. "r") is flattend.
% The belief about the volatility of the outcome and the volatility of the
% volatility is not altered (as these aren't necessarily linked to specific
% stimuli). Input arguments are "info" which is a vector (of length ntrials) of 1s or 0s coding
% for the binary outcome being learned about."flattenpoints" is a vector of
% trial numbers at which the flattening occurs. "out" is a structure
% containing the relevant output of hte model (beliefs and distributions
% across all three levels, that sort of thing).

if nargin<2
    flattenpoints=0;
end
if ~isvector(info) || ~isnumeric(info) || min(info)<0 || max(info)>1
    error('Info should be a vector of 0s or 1s')
end

if ~isvector(flattenpoints) || ~isnumeric(flattenpoints) || min(flattenpoints)<0 || max(flattenpoints)>length(info)
   error('Flattenpoints should be a vector of trial numbers with each number > 1 and < total number of trials in info')
end



% r-axis (reward rate)
rvec = .01:.02:.99;
rSize = length(rvec);

% Ilog axis (I = 1/v)
%Note minimum precision (1/v) is 2 as is defined as a+b from the beta
%distribution and the initial prior is a=1, b=1 (i.e. flat)
Ilog = log(2):0.2:log(10000);
Isize = length(Ilog);

% klog axis
klog = log(5e-4):.2:log(2000);
kSize = length(klog);

tmp = ones(rSize, Isize)*(1/(rSize*Isize*kSize));
% load tmp into rIk for each k
% pIk: 3D distribution that we're always updating p(r_i, I_i, k | every
% observation to date).
rIk = zeros(rSize, Isize, kSize);
for k = 1 : kSize
    rIk(:,:,k) = tmp;
end
rIk_old=rIk;

% precompute p(Ii+1|Ii,k) for every I_{i+1}, I_i, k
Ip1gIk = zeros(Isize, Isize, kSize); % 3d array for I transition.
tmpI = zeros(Isize, Isize);  % to keep N(I_{i+1},k2
tmpI2=tmpI;
for k = 1 : kSize
    for I = 1:Isize
        for Ip1 = 1:Isize
            % N(I_{i+1},k2)
            var = exp(klog(k)*2); % k is stdev
            % below is the probability density function of a gaussian with
            % mean Ip1 and SD k at point I. It essentially provides the probability of
            % Ip1 given I and k
          % tmpI(Ip1, I) = (exp(-power((Ilog(I) - Ilog(Ip1)),2))/(2*var)) / (sqrt(2*pi*var)); % there is an error in the brackets here (the (2*var) should be part of the exponential) 
            tmpI(Ip1, I) = (exp(-power((Ilog(I) - Ilog(Ip1)),2)/(2*var))) / (sqrt(2*pi*var)); % this is corrected
        end
        % normalise so p(i_i+1|I_i,k) sums to 1; 
        tmpI(:,I) = tmpI(:,I)./sum(tmpI(:,I)); % this is a right stochastic (transitional) matrix. 
    end
   
    %Probability of I on next trial given I and k on current trial. 
    Ip1gIk(:,:,k) = tmpI; % place tmpI in it. (overall this should sum to ksize*Isize). 
end


% precompute p(ri+1|ri,Ii+1) for every r_{i+1}, r_i, I_{i+1}
rp1gpIp1 = zeros(rSize, rSize, Isize);
tmpp = zeros(rSize, rSize);
for Ip1 = 1:Isize
    for r = 1:rSize
        for rp1 = 1:rSize
            % Beta(x; a,b) NB can think of this where a is number of heads
            % and b number of tails. The initial setting is a=1 and b=1 which
            % gives a flat prior, the mean and SD add to this
            %
            a = 1 + (exp(Ilog(Ip1)))*rvec(r);
            b = 1 + (exp(Ilog(Ip1)))*(1-rvec(r));
            x = rvec(rp1);
            if ~(x==0 || x==1)  % changed from ~(x==0)||~(x==1), which is wrong. beta pdf not defined for 0 and 1
                %NB beta pdf is ((x.^(a-1)).*((1-x).^(b-1)))./denominator
                %(betaln_ab below). This code works out log of this
                logkerna = (a-1)*log(x);
                logkernb = (b-1)*log(1-x);
                betaln_ab = gammaln(a) + gammaln(b) - gammaln(a+b); % log denomintar of the beta pdf
                tmpp(rp1, r) = exp(logkerna + logkernb - betaln_ab); %beta pdf (no longer in log space)
            else
                tmpp(rp1, r) = 0;
            end
            
        end
        tmpp(:,r) = tmpp(:,r)./sum(tmpp(:,r)); % normalise across each range of pi+1
    end
    rp1gpIp1(:,:,Ip1) = tmpp; % place tmpp in it. Probability of r on next trial given r on this trial and I on next trial
end

%
% BAYESIAN UPDATE
%
% this is liklihood across all trials in task (without taking into account the volatility or k). Not
% really sure why the k loop is needed
for trial = 1:length(info)
    
      if sum(trial==flattenpoints)>0
      rIk=repmat(mean(rIk),rSize,1,1);
        rIk_old=rIk; 
      end
    
      
      %entropy
      out.entropy(trial,1)=sum(sum(rIk(rIk>eps).*(log2(1./(rIk(rIk>eps))))));  % numbers less than 1e-15 are 0
    
        % GET MARGINALS
    %
    
    % p(reward)
    pI = sum(rIk,3);
    out.pDist(trial,:) = sum(pI,2); % sum over each row.
    out.pEst(trial,:)  = sum(out.pDist(trial,:).*rvec);
    
    % p(volatility)
    pI = sum(rIk,3);
    out.iDist(trial,:) = sum(pI); % column-wise sum.
    out.iEst(trial,:)  = sum(out.iDist(trial,:).*Ilog); % note this returns preciscion
    out.vEst(trial,:)=1./(out.iEst(trial,:)); % this is the actual volatility

    
    %k(volatility of volatility)
    out.kDist(trial,:)=sum(sum(rIk,1),2);
    out.kEst(trial,:)=sum(out.kDist(trial,:).*klog);
    
    
    
    
    if (info(trial)==1)
       % for k = 1:kSize
            for r=1:rSize
                rIk(r,:,:) = rIk(r,:,:)*rvec(r);
            end
        %end
    else
        %for k = 1:kSize
            for r=1:rSize
                rIk(r,:,:) = rIk(r,:,:)*(1-rvec(r));
            end
        %end
        
    end
    
    % now do normalization
    rIk = rIk ./ sum(sum(sum(rIk)));
    
    %
    % INFORMATION LEAK (increase the variance in the joint distribution)
    %
    
    % I) multiply rIk (simple likelihood) by Ip1gIk (probability of I on 
    %next trial given I and k on this trial), and integrate out I. This 
    %will give pIp1k (probability of r on next trial given I+1).
    for k = 1:kSize
        pIp1k = zeros(rSize,Isize);
        for Ip1 = 1:Isize
            for r = 1:rSize
                pIp1k(r,Ip1) = sum(Ip1gIk(Ip1,:,k).*rIk(r,:,k)); % for a given k calcuate the probability of r given the updated I
            end
        end
        % II) multiply pIp1k (probability of r on next trial given I) by 
        %pp1gpIp1 (probability of r on next trial given r on this trial and 
        %I on next trial), and integrate out r on this trial. This will give pp1Ip1k 
        %(probability of r on next trial given I on nextr trial and k).
        pp1Ip1k = zeros(rSize,Isize);
        for Ip1 = 1:Isize
            for rp1 = 1:rSize
                pp1Ip1k(rp1,Ip1) = sum(pIp1k(:,Ip1).*rp1gpIp1(rp1,:,Ip1)');
            end
        end
        % III) Place pp1Ip1k into pIk (belief that is carried to the next
        % trial.
        rIk(:,:,k) = pp1Ip1k;
    end
     out.KLdiv(trial,1)=sum(sum(rIk(rIk>eps & rIk_old>eps).*(log2(rIk(rIk>eps & rIk_old>eps))-log2(rIk_old(rIk>eps & rIk_old>eps)))));   % definition of kldiv when prob is 0 is 0-- give identical value to above (but is faster).
   rIk_old=rIk;
    %

    
    
 end

out.rvec=rvec;
out.Ilog=Ilog;
out.klog=klog;
out.info=info;
out.flattenpoints=flattenpoints;






