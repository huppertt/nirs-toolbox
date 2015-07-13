classdef FixSatChans < nirs.modules.AbstractModule

    properties
    end
    
    methods

        function obj = FixSatChans( prevJob )
           obj.name = 'Fix Oversaturated Channels';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
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
                    a(j,1) = length(ar_fit(d(:,lst(j)), round(10*Fs)));
                end
                
                % bad chans = high variance and autcorr
                bad = lst(a > 3*Fs);
                
                d(:,bad) = lognrnd( 0, 1, [size(d,1) length(bad)] );
                
                data(i).data = d;
            end
        end
    end
    
end

