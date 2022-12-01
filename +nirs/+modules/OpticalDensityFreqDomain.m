classdef OpticalDensityFreqDomain < nirs.modules.AbstractModule
%% OpticalDensity - Converts raw data to optical density using FD-NIRS data.
% 
% dOD = -log( raw/mean(raw) )

methods
    function obj = OpticalDensityFreqDomain ( prevJob )
        obj.name = 'Optical Density Frequency Domain ';
        if nargin > 0
            obj.prevJob = prevJob;
        end
    end
    
    function data = runThis( obj, data )
        for i = 1:numel(data)
            if(~ismember('MeasurementGroup',data(i).probe.link.Properties.VariableNames))
                warning(['Measurement groups not defined- skipping file']);
                continue;
            end
            newlink=table;
            newd=[];
            Groups=data(i).probe.link.MeasurementGroup;
            uGroups=unique(Groups);
            for gIdx=1:length(uGroups)
                lst=find(Groups==uGroups(gIdx));
                utypes=unique(data(i).probe.link(lst,:).type);
                for tIdx=1:length(utypes)
                    
                    lst2=find(Groups==uGroups(gIdx) & data(i).probe.link.type==utypes(tIdx));
                    modfreq=unique(data(i).probe.link.ModFreq(lst2));
                    for mIdx=1:length(modfreq)
                        lst2=find(Groups==uGroups(gIdx) & data(i).probe.link.type==utypes(tIdx) & ...
                            data(i).probe.link.ModFreq==modfreq(mIdx));
                        d=data(i).data(:,lst2);
                        dist=data(i).probe.distances(lst2);
                        
                        [absorp,scatter]=computeFDmultidist(d,dist,modfreq(mIdx));
                        ROI=uGroups(gIdx);
                        if(true)
                            newd=[newd; absorp'];
                            newlink=[newlink; table(ROI) data(i).probe.link(lst2(1),:)];
                        else
                            newd=[newd; scatter'];
                            newlink=[newlink; table(ROI) data(i).probe.link(lst2(1),:)];
                        end
                        
                    end
                end
            end
            newlink.source=[];
            newlink.detector=[];
            newlink.ModFreq=[];
            newlink.MeasurementGroup=[];
            newlink.Calibration_Applied=[];
            newlink.Calibration_Factor=[];
            newlink.Calibration_Term=[];
            
            probeROI=nirs.core.ProbeROI;
            probeROI.link=newlink;
            probeROI.RegionNames=cellstr(num2str(uGroups));
            data(i).probe=probeROI;
            data(i).data=newd';
        end
    end
    
end
end


function [absoprtion,scattering]=computeFDmultidist(d,dist,ModFreq)
wmod = 2*pi*ModFreq*1E6;
v = 2.26100e+10; 

ac = abs(d);
phs = angle(d);
X(:,1)=dist;
X(:,2)=1;


Uac = log(ac.*(ones(size(ac,1),1)*dist'.^2));
u=Uac./(ones(size(Uac,1),1)*median(Uac,1));
w = inv(chol(cov(u)));
wX=w*X;
Sac=inv(wX'*wX)*wX'*w*Uac';
Sac =Sac(1,:)';

p=phs./(ones(size(phs,1),1)*median(phs,1));
w = inv(chol(cov(p)));
wX=w*X; %(:,1);
Sph = inv(wX'*wX)*wX'*w*phs';
Sph =Sph(1,:)';


absoprtion= wmod/(2*v)*(Sph./Sac-Sac./Sph);   % Eqn 10 from Fantini 1994
scattering= (Sac.^2-Sph.^2)./(3*absoprtion)-absoprtion;

end
