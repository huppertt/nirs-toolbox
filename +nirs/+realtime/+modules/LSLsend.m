classdef LSLsend < nirs.realtime.modules.AbstractModule
    % real-time implementation to send data across LSL
    
    properties
      dataExtractionFcn; % this function needs to have the form val = fcn(d,t,probe,stimulus)
      LSLdata_StreamName='NIRStoolbox';
      LSLdatatype='data';  % or marker
      LSLdataNumChannels=1;
    end
    properties(Hidden=true)
        LSLdata_Stream=[];
        liblsl=[];
      
    end
    
    methods
        function obj=LSLsend(prevJob)
            obj.name='RT-LSL sender';
            if nargin > 0
                obj.prevJob = prevJob;
            end
              obj.liblsl = lsl_loadlib();
              
             
              
        end
        
        function obj=resetThis(obj)
             obj.LSLdata_Stream= [];
              
        end
        
        
        function [dOut,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(isempty(obj.LSLdata_Stream))
                 if(strcmp(obj.LSLdatatype,'data'))
                info = lsl_streaminfo(obj.liblsl,obj.LSLdata_StreamName,'EEG',obj.LSLdataNumChannels,0,'cf_float32','nirsrealtimemodulesLSLsend');
              else
                      info = lsl_streaminfo(obj.liblsl,obj.LSLdata_StreamName,'Markers',obj.LSLdataNumChannels,0,'cf_string','nirsrealtimemodulesLSLsend');
              end
                obj.LSLdata_Stream= lsl_outlet(info);
            end
            
            dOut=d;
            val=feval(obj.dataExtractionFcn,d,t,probe,stimulus);
            
            if(strcmp(obj.LSLdatatype,'data'))
               obj.LSLdata_Stream.push_sample(val);
            else
                if(isnumeric(val))
                        val={num2str(val)};
                end
                obj.LSLdata_Stream.push_sample({val});
            end
            
            
        end
            
    end
    
end