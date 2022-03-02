classdef FixFlatChans < nirs.modules.AbstractModule
%% FixFlatChans - Replaces flatlined channels with high variance noise.
%
% Note: This should be run on raw data (before optical density conversion).
    properties
        window_length=3; % smoothing window
        remove_full_channel=true;
    end
    methods
        function obj = FixFlatChans( prevJob )
           obj.name = 'Fix Flatlined Channels';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                
                d = data(i).data;
                Fs = data(i).Fs;
                
                % subtract mean
                m= mean(d,1);
                v = bsxfun(@minus, d,m);
                
                v = v.^2;
                
                % 3s windowed variance
                f = round(Fs*obj.window_length);
                v = filter(f, 1, v);
                
                
                
                for id=1:size(d,2)
                    % fill with noise
                    if(obj.remove_full_channel)
                        % bad channels
                        if(find(any(v(:,id) < 1e-9)))
                            d(:,id) = NaN;
                        end
                    else
                        lst=find(v(:,id) < 1e-9);
                        d(lst,id)=NaN;
                        
                    end
                end
                data(i).data = d;
            end
        end
    end
    
end

