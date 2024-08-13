function rts = cal_RTs(RT,winchosen,losschosen,choice)
% this sorts and medians RTs for different different trials types

%median RT for the block
rts.all=median(RT);
%median RT for trials where win received or not received in the previous trials
winlasttrial=nan(size(winchosen));
winlasttrial(2:end)=winchosen(1:end-1);
rts.afterwin=median(RT(winlasttrial==1));
rts.afternowin=median(RT(winlasttrial==0));
%median RT for trials where loss received or not received in the previous trials
losslasttrial=nan(size(losschosen));
losslasttrial(2:end)=losschosen(1:end-1);
rts.afterloss=median(RT(losslasttrial==1));
rts.afternoloss=median(RT(losslasttrial==0));
%median RT for trials where a switch has been made or not been made
ifswtich=nan(size(choice));
ifswtich(2:end)=choice(2:end)-choice(1:end-1);
rts.stay=median(RT(ifswtich==0));
rts.switch=median(RT(abs(ifswtich)==1));
end

