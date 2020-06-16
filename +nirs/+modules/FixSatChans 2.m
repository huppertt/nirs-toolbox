classdef FixSatChans < nirs.modules.AbstractModule
%% FixSatChans - Attempts to replace oversaturated channels with high variance noise.
%
% Note: This should be run on raw data (before optical density conversion).
    
    methods
        function obj = FixSatChans( prevJob )
           obj.name = 'Fix Oversaturated Channels';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                
                d = data(i).data;
                Fs = data(i).Fs;
                
                % variance
                v = std(d,1)';
                
                u = median(v);
                s = mad(v,0) / 0.6745;
                
                t = (v-u)/s;
                
                lst = find( abs(t) > 2 );
                
                a = zeros( size(lst) );
                for j = 1:length(lst)
                    a(j,1) = length(nirs.math.ar_fit(d(:,lst(j)), round(10*Fs)));
                end
                
                % bad chans = high variance and autocorr
                bad = lst(a > 3*Fs);
                if(~isempty(bad))
                    m= mean(d,1);
                    d(:,bad) = lognrnd( 0, 1, [size(d,1) length(bad)] );
                    d(:,bad) = bsxfun(@plus, d(:,bad),m(bad));
                    data(i).data = d;
                end
            end
        end
    end
    
end

