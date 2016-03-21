classdef CalculateCMRO2 < nirs.modules.AbstractModule
    %% CalculateCMRO2 - This function will add CMRO2 and Flow
    %  The input must have the hbo/hbr data types.
    %
    
    properties
        HbT_baseline=50; % baseline uM HbT
        E0=.4;
        tau=3;
        alpha=.4;
        
      
    end
    
    methods
        
        function obj = CalculateCMRO2( prevJob )
            obj.name = 'Add CMRO2 and CBF to data using a Kalman filter';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            if(~all(ismember({'hbo','hbr'},data(1).probe.link.type)))
                warning('data does not contain oxy/deoxy-Hb.  Use the MBLL first');
                return
            end
            
            nlgr =WKM_idnl(obj.HbT_baseline,obj.E0,obj.tau,obj.alpha);
            
            for i = 1:length(data)
                %Time series models
                link=data(i).probe.link;
                link=link(ismember(link.type,'hbo'),:);
                CMRO2=zeros(size(data(i).data,1),height(link));
                CBF=zeros(size(data(i).data,1),height(link));
                
                for j=1:height(link)
                    iHbO=find(ismember(data(i).probe.link(:,1:2),link(j,1:2)) & ismember(data(i).probe.link.type,'hbo'));
                    iHbR=find(ismember(data(i).probe.link(:,1:2),link(j,1:2)) & ismember(data(i).probe.link.type,'hbr'));
                    HbO2=data(i).data(:,iHbO);
                    HbR=data(i).data(:,iHbR);
                    
                    d=iddata([HbO2,HbR],[],1./data(i).Fs,'OutPutName',{'HbO2','HbR'});
                    
                    [CMRO2(:,j),CBF(:,j)]=kalman_fit(nlgr,d);
                    CMRO2(:,j)=CMRO2(:,j)*100;% convert to % signal change so the scale is closer to the HB signals
                    CBF(:,j)=CBF(:,j)*100;
                end
                linkCMRO2=link;
                linkCMRO2.type=repmat({'CMRO2'},height(link),1);
                linkCBF=link;
                linkCBF.type=repmat({'CBF'},height(link),1);
                data(i).data=[data(i).data CMRO2 CBF];
                data(i).probe.link=[data(i).probe.link; linkCBF; linkCMRO2];
                data(i)=sorted(data(i));
                
                
            end
        end
    end
    
end

