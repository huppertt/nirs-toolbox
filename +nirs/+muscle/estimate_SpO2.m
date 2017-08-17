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

for i=1:height(link)
    chaIdx=find(b==i);
    R=(ub(:,chaIdx)-lb(:,chaIdx))./mean(.5*(ub(:,chaIdx)+lb(:,chaIdx)));
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
        
    SpO2(:,i)=(E(1,2)-R*E(2,2))./(R*(E(2,1)-E(2,2))+(E(1,2)-E(1,1)));
end

SpO2=min(SpO2,1);
SpO2=max(SpO2,0);