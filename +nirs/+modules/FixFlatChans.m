classdef FixFlatChans < nirs.modules.AbstractModule
%% FixFlatChans - Replaces flatlined channels with high variance noise.
%
% Note: This should be run on raw data (before optical density conversion).
    
    methods
        function obj = FixFlatChans( prevJob )
           obj.name = 'Fix Flatlined Channels';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                Fs = data(i).Fs;
                
                % subtract mean
                v = bsxfun(@minus, d, mean(d,1));
                
                v = v.^2;
                
                % 3s windowed variance
                f = round(Fs*3);
                v = filter(f, 1, v);
                
                % bad channels
                bad = find(any(v < 1e-9));
                
                % fill with noise
                d(:,bad) = lognrnd(0, 1, [size(d,1) length(bad)]);
                data(i).data = d;
            end
        end
    end
    
end

