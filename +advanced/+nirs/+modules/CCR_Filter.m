classdef CCR_Filter < nirs.modules.AbstractModule
%% Multual information filter based on canonical correlation of 

    
    properties
       cutoff = 0.05;   % statistical cutoff
    end
    
    methods

        function obj = CCR_Filter( prevJob )
           obj.name = 'CCR Filter';
           if nargin > 0
               obj.prevJob = prevJob;
           end
            
        end
        
        function data = runThis( obj, data )
            for idx = 1:numel(data)
                d = data(idx).data;
   
                types=unique(data(idx).probe.link.type);
                if(~iscell(types))
                    types=num2cell(types);
                end
                
                if(length(types)~=2)
                    error('This module only supports two data types currently');
                end
                
                D1=data(idx).data(:,ismember(data(idx).probe.link.type,types{1}));
                D2=data(idx).data(:,ismember(data(idx).probe.link.type,types{2}));
                
                % [D1,f1] = nirs.math.innovations(D1,data(idx).Fs*4,true);
                % [D2,f2] = nirs.math.innovations(D2,data(idx).Fs*4,true);
           
                 [A,B,R,U,V,stats]=canoncorr(D1,D2);
                 
                 lst=find(nirs.math.fdr(stats.p)<0.05);
                 
%                  if length(lst) < 4
%                      lst = [1 2 3];
%                  end
%                  
                D1 = U(:,lst)*pinv(A(:,lst))+ones(size(D1,1),1)*mean(D1,1);
                D2 = V(:,lst)*pinv(B(:,lst))+ones(size(D2,1),1)*mean(D2,1);
                
                
%                 for i=1:length(f1)
%                     D1(:,i)=filter(1,[1; f1{i}(2:end)],D1(:,i));
%                 end
%                 for i=1:length(f2)
%                     D2(:,i)=filter(1,[1; f2{i}(2:end)],D2(:,i));
%                 end
                
                data(idx).data(:,ismember(data(idx).probe.link.type,types{1}))=D1;
                data(idx).data(:,ismember(data(idx).probe.link.type,types{2}))=D2;
                 
            end
        end
    end
    
end

