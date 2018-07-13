classdef CalculateTotalHb < nirs.modules.AbstractModule
    %% CalculateTotalHb - This function will add total-hemoglobin and
    % the tissue oxygen index estimates to the data structure.  This function
    % will accept either nirs.core.Data (time-series) or nirs.core.ChannelStats
    % variables.  The input must have the hbo/hbr data types. 
    %
    % NOTE- Tissue oxygen saturation is defined as:
    %  [HbO2(t=0)+dHbO2(t)]/[HbO2(t=0)+HbR(t=0)+dHbO2(t)+dHbR(t)]
    
    properties
        StO2_baseline=.70;
        HbT_baseline=50; % baseline uM HbT
        PPF=.1;  % PPF term used in the MBLL (because I need to remove this for the SO2 calculation)
    end
    
    methods
        
        function obj = CalculateTotalHb( prevJob )
            obj.name = 'Add total Hb and TOI to data';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            if(~all(ismember({'hbo','hbr'},data(1).probe.link.type)))
                warning('data does not contain oxy/deoxy-Hb.  Use the MBLL first');
                return
            end
            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.Data'))
                    %Time series models
                    link=data(i).probe.link;
                    link=link(ismember(link.type,'hbo'),:);
                    HbT=zeros(size(data(i).data,1),height(link));
                    TOI=zeros(size(data(i).data,1),height(link));
                    
                    for j=1:height(link)
                        iHbO=find(ismember(data(i).probe.link(:,1:2),link(j,1:2)) & ismember(data(i).probe.link.type,'hbo'));
                        iHbR=find(ismember(data(i).probe.link(:,1:2),link(j,1:2)) & ismember(data(i).probe.link.type,'hbr'));
                        HbO2=data(i).data(:,iHbO);
                        HbR=data(i).data(:,iHbR);
                        HbT(:,j)=HbO2+HbR;
                        TOI(:,j)=(obj.HbT_baseline.*obj.StO2_baseline+HbO2*obj.PPF)./max(obj.HbT_baseline+HbO2*obj.PPF+HbR*obj.PPF,1);
                        
                    end
                    linkTOI=link;
                    linkTOI.type=repmat({'StO2'},height(link),1);
                    linkHbT=link;
                    linkHbT.type=repmat({'hbt'},height(link),1);
                    data(i).data=[data(i).data HbT TOI];
                    data(i).probe.link=[data(i).probe.link; linkHbT; linkTOI];
                    data(i)=sorted(data(i));
                elseif(isa(data(i),'nirs.core.ChannelStats'))
                    variab=data(i).variables;
                    variab=variab(ismember(variab.type,'hbo'),:);
                    lst=ismember(variab.Properties.VariableNames,{'source','detector','cond'});
                    for j=1:height(variab)
                        iHbO=find(ismember(data(i).variables(:,lst),variab(j,lst)) & ismember(data(i).variables.type,'hbo'));
                        iHbR=find(ismember(data(i).variables(:,lst),variab(j,lst)) & ismember(data(i).variables.type,'hbr'));
                        
                        %HbT is easy
                        c=zeros(1,height(data(i).variables));
                        c([iHbO iHbR])=1;
                        betaHbT(j,1)=c*data(i).beta*obj.PPF;
                        CovHbT(j,j)=c*data(i).covb*c'*obj.PPF^2;
                        
                        % SO2 is much harder
                        %SO2=(obj.HbT_baseline.*obj.StO2_baseline+HbO2*obj.PPF)./max(obj.HbT_baseline+HbO2*obj.PPF+HbR*obj.PPF,1)-obj.StO2_baseline;
                        betaSO2(j,1)=(obj.HbT_baseline.*obj.StO2_baseline+data(i).beta(iHbO)*obj.PPF)./...
                            max(obj.HbT_baseline+data(i).beta(iHbO)*obj.PPF+data(i).beta(iHbR)*obj.PPF,1)-obj.StO2_baseline;
                        
                        % (HbO0+dHbO)/(HbO0+dHbO+HbR0+dHbR) =>
                        % [(HbO0+dHbO+HbR0+dHbR)*SE(HbO2) -
                        % (HbO0+dHbO)*SE(HbT)]/ (HbO0+dHbO+HbR0+dHbR)^2
                        SEHbO=sqrt(data(i).covb(iHbO,iHbO))*obj.PPF;
                        SEHbT=sqrt(CovHbT(j,j));
                        HbO=(obj.HbT_baseline.*obj.StO2_baseline+data(i).beta(iHbO)*obj.PPF);
                        HbT=max(obj.HbT_baseline+data(i).beta(iHbO)*obj.PPF+data(i).beta(iHbR)*obj.PPF,1);
                        
                        CovSO2(j,j)= ((HbT*SEHbO+HbO*SEHbT)./HbT.^2).^2;
                        
                    end
                    variabTOI=variab;
                    variabTOI.type=repmat({'StO2'},height(variab),1);
                    variabHbT=variab;
                    variabHbT.type=repmat({'hbt'},height(variab),1);
                    
                    data(i).variables=[data(i).variables; variabHbT; variabTOI];
                    data(i).beta=[data(i).beta; betaHbT; betaSO2];
                    data(i).covb=blkdiag(data(i).covb,CovHbT,CovSO2);  % This isn't entirely correct, but ok for now
                    
                    link=data(i).probe.link;
                    link=link(ismember(link.type,'hbo'),:);
                    linkTOI=link;
                    linkTOI.type=repmat({'StO2'},height(link),1);
                    linkHbT=link;
                    linkHbT.type=repmat({'hbt'},height(link),1);
                    data(i).probe.link=[data(i).probe.link; linkHbT; linkTOI];
                    
                    data(i)=sorted(data(i));
                    
                end
                
            end
        end
    end
    
end

