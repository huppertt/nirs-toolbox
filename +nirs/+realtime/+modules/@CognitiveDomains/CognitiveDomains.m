classdef CognitiveDomains < nirs.realtime.modules.AbstractModule
    properties
        model={'mindfulness.nii.gz'}
        training=120;  % number of samples to delay the model to get an estimate of the noise
        LSLdata_StreamName='NIRStoolbox';
        addplot=true;
        
    end
    properties(Hidden=true)
        projector=[];
        kalmanfilter=[];
        trainingdata=[];
        iVar=[];
        LSLdata_Stream=[];
        liblsl=[];
        baseline_val=[];
        linehandles=[];
    end
    methods
        function obj = CognitiveDomains(prevJob)
            obj.name='RT cognitive domains model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
            obj.liblsl = lsl_loadlib();
        end
        
        function obj=resetThis(obj)
            obj.trainingdata=[];
            obj.projector=[];
            obj.LSLdata_Stream= [];
            obj.iVar=[];
            obj.baseline_val=[];
        end
        
        function [d,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(isempty(obj.projector))
                disp('Gathering training data for CognitiveDomains model');
                for i=1:length(obj.model)
                    p=fileparts(which('nirs.realtime.modules.CognitiveDomains'));
                    roi=nirs.util.convertNIFTI_to_roi(probe,fullfile(p,obj.model{i}));
                    obj.projector(i,:)=roi.weight(1:2:end);
                end
                obj.kalmanfilter=nirs.realtime.util.RobustKalmanFilter(1E-14);
                info = lsl_streaminfo(obj.liblsl,obj.LSLdata_StreamName,'EEG',length(obj.model),0,'cf_float32','nirsrealtimemodulesLSLsend');
                obj.LSLdata_Stream= lsl_outlet(info);
            end
            
            if(size(obj.trainingdata,1)<obj.training)
                obj.trainingdata = [obj.trainingdata; d];
            else
                if(isempty(obj.iVar))
                    obj.iVar=sqrt(inv(diag(var(obj.trainingdata(:,1:2:end),[],1))));
                    obj.projector=obj.projector*obj.iVar;
                    for i=1:size(obj.trainingdata,1)
                        obj.kalmanfilter.update(abs(obj.projector*obj.trainingdata(i,1:2:end)'),1,1);
                    end
                    val=zeros(length(obj.model),1);
                    for i=1:size(obj.trainingdata,1)
                        obj.kalmanfilter.update(abs(obj.projector*obj.trainingdata(i,1:2:end)'),1,1);
                        val=val+obj.kalmanfilter.B./sqrt(diag(obj.kalmanfilter.P));  % ttest value
                    end
                    obj.baseline_val=val/size(obj.trainingdata,1);
                    
                    if(obj.addplot)
                        figure; hold on;
                        for i=1:length(val)
                            obj.linehandles(i)=plot(NaN,NaN);
                        end
                        legend(obj.model)
                    end
                    
                    disp('Finished training CognitiveDomains model');
                else
                    
                    for i=1:size(d,1)
                        obj.kalmanfilter.update(abs(obj.projector*d(i,1:2:end)'),1,1);
                        val=obj.kalmanfilter.B./sqrt(diag(obj.kalmanfilter.P))-obj.baseline_val;  % ttest value
                        obj.LSLdata_Stream.push_sample(val);
                        
                        
                        if(obj.addplot)
                            for j=1:length(val)
                                xdata=[get(obj.linehandles(j),'Xdata') t(i)];
                                ydata=[get(obj.linehandles(j),'Ydata') val(j)];
                                set(obj.linehandles(i),'Xdata',xdata(:),'Ydata',ydata(:));
                            end

                        end
                        
                    end
                end
            end
            
        end
        
    end
    
    
end