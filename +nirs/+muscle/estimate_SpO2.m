function SpO2 = estimate_SpO2(rawInten)
%
% Calibration-Free Pulse Oximetry Based on Two Wavelengths in
% the Infrared ? A Preliminary Study
% Meir Nitzan, Salman Noach, Elias Tobal, Yair Adar, Yaacov Miller, Eran Shalom 
% and Shlomo Engelberg. Sensors 2014, 14, 7420-7434; doi:10.3390/s140407420


if(length(rawInten)>1)
    for i=1:length(rawInten)
        SpO2(i)=nirs.muscle.estimate_SpO2(rawInten(i));
    end
    return
end


link=rawInten.probe.link;
link.type=[];
[link,~,b]=unique(link);

[lb,ub]=nirs.muscle.findPulseAC(rawInten.data,rawInten.Fs);

if(iscell(rawInten.probe.link.type))
for i=1:height(rawInten.probe.link)
    type(i,1)=str2num(rawInten.probe.link.type{i});
end
rawInten.probe.link.type=type;
end

for i=1:height(link)
    chaIdx=find(b==i);
    R=(ub(:,chaIdx)-lb(:,chaIdx))./(ones(size(ub,1),1)*mean(.5*(ub(:,chaIdx)+lb(:,chaIdx)),1));
    if(size(R,2)~=2)
        error('not yet supported')
    else
        R=R(:,1)./R(:,2);
    end
    
    lambda = rawInten.probe.link.type(chaIdx);
    ext = nirs.media.getspectra( lambda );
    
    clist = [1 2]; % hbo and hbr; need to fix this
    
    % extinction coefficients
    E = ext(:,clist);
        
    SO2(:,i)=(E(1,2)-R*E(2,2))./(R*(E(2,1)-E(2,2))+(E(1,2)-E(1,1)));
end

SO2=min(SO2,1);
SO2=max(SO2,0);

SpO2 = rawInten;
SpO2.data=SO2;
SpO2.probe.link=link;
SpO2.probe.link.type=repmat({'SpO2'},height(link),1);




