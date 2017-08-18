function pFlow = estimate_pFlow(hb)
%
% Themelis, G., D?Arceuil, H., Diamond, S. G., Thaker, S., Huppert, T. J., Boas, D., & Franceschini, M. A. 
% (2006, March). NIRS Measurement of the Pulsatile Component of Cerebral Blood Flow 
% and Volume from the Arterial Oscillations. In Biomedical Topical Meeting (p. SH26). Optical Society of America.

if(length(hb)>1)
    for i=1:length(hb)
        pCBF(i) = nirs.muscle.estimate_rFlow(hb(i));
    end
    return;
end

lst=find(ismember(hb.probe.link.type,'hbo'));
lst2=find(ismember(hb.probe.link.type,'hbr'));

if(isempty(lst))
    error('data must contain HbO2 and Hb variables');
end

d=hb.data(:,lst)+hb.data(:,lst2);

[fa,fb]=butter(4,[.2]*2/hb.Fs,'high');
d=filtfilt(fa,fb,d);


%find the upper/lower bound of the cardiac cycle
[lb,ub]=nirs.muscle.findPulseAC(d,hb.Fs);

%midpoint
mid=(lb+ub)/2;

for i=1:size(d,2)

    %find the time-points when the cardiac crosses the mid-point 
    lstCross = find(d(1:end-1,i)<mid(1:end-1,i) & d(2:end,i)>=mid(2:end,i));
    lstCross=[1; lstCross; size(d,1)];
    rF=[];
    t=[];
    for j=2:length(lstCross)-1
        lstN=lstCross(j-1):lstCross(j);
        lstP=lstCross(j):lstCross(j+1);

        %for each crossing, find the highest and lowest point 
        [n,iN]=min(d(lstN,i));
        [p,iP]=max(d(lstP,i));
        rF(j)=(p-n)/(hb.time(lstP(iP))-hb.time(lstN(iN)));
        
        t(j)=hb.time(lstCross(j));
    end
    CBF(:,i)=interp1(t(~isnan(rF)),rF(~isnan(rF)),hb.time,'spline','extrap');
    CBF(:,i)=medfilt1(CBF(:,i),fix(6*hb.Fs));
end

CBF=CBF./(ones(size(CBF,1),1)*median(CBF,1));

pFlow = hb;
pFlow.data=CBF;
pFlow.probe.link=hb.probe.link(lst,:);
pFlow.probe.link.type=repmat({'pFlow'},length(lst),1);





