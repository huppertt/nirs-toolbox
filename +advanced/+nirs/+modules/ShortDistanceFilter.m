classdef ShortDistanceFilter < nirs.modules.AbstractModule
%% ShortDistanceFilter - This function filters based on short distance measurements
%
% Options:
%     ncomp - % number of components to remove

    properties
         maxnumcomp=6;
         splittypes=false;
    end
    
    methods

        function obj =  ShortDistanceFilter( prevJob )
           obj.name = 'Short Distance Correction';         
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
           
           
            for i = 1:numel(data)
                disp([num2str(i) ' of ' num2str(numel(data))])
                if(obj.splittypes)
                    types=unique(data(i).probe.link.type);
                    for tI=1:length(types)
                        channelLst=find(data(i).probe.link.ShortSeperation & ismember(data(i).probe.link.type,types(tI)));
                        lst=find(~data(i).probe.link.ShortSeperation & ismember(data(i).probe.link.type,types(tI)));
                        d = data(i).data;
                        
                        SS = d(:,channelLst);
                        if(obj.maxnumcomp<=size(SS,2))
                            [U,S,V]=nirs.math.mysvd(SS);
                            SS=U(:,1:obj.maxnumcomp)*S(1:obj.maxnumcomp,1:obj.maxnumcomp)*V(:,1:obj.maxnumcomp)';
                        end
                        SS=orth(SS);
                        SS(:,end+1)=1;
                        
                        b = nirs.math.ar_irls(d(:,lst),SS,fix(data(i).Fs*4));
                        data(i).data(:,lst)=data(i).data(:,lst)-SS*b.beta;
                    end
                    
                else
                    channelLst=find(data(i).probe.link.ShortSeperation);
                    lst=find(~data(i).probe.link.ShortSeperation);
                    d = data(i).data;
                    
                    SS = d(:,channelLst);
                    if(obj.maxnumcomp<=size(SS,2))
                        [U,S,V]=nirs.math.mysvd(SS);
                        SS=U(:,1:obj.maxnumcomp)*S(1:obj.maxnumcomp,1:obj.maxnumcomp)*V(:,1:obj.maxnumcomp)';
                    end
                    SS=orth(SS);
                    SS(:,end+1)=1;
                    
                    b = nirs.math.ar_irls(d(:,lst),SS,fix(data(i).Fs*4));
                    data(i).data(:,lst)=data(i).data(:,lst)-SS*b.beta;
                end
               
                
            end
        end
    end
    
end

