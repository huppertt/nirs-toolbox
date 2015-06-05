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
                
                p = tcdf(-abs(v-u)/s, length(v)-1);
                
%                 % power spec
%                 a = fft(d);
%                 a = abs( a(2:end,:) );
%                 
%                 u = median(a)';
%                 s = mad(a,0)' / 0.6745;
%                 
%                 t = bsxfun(@minus, a, u');
%                 t = bsxfun(@rdivide, t, s');
%                 
%                 p2 = tcdf(-abs(t), length(a)-1);
%                 p2 = sum( p2 < 1e-6, 1 )';

%                 
%                 % autocorr
%                 a = zeros( size(d,2), round(4*Fs)+1 );
%                 for j = 1:size(d,2)
%                     a(j,:) = autocorr(d(:,j), round(4*Fs));
%                 end
%                 
%                 a = sum(abs(a),2);

                for j = 1:size(d,2)
                    a(j,1) = length(ar_fit(d(:,j), round(10*Fs)));
                end
                
                % bad chans = high variance and autcorr
                lst = (p < 0.05) & (a > 3*Fs);
                
                d(:,lst) = lognrnd( log(100), log(1.1), [size(d,1) sum(lst)] );
                
                disp([i sum(lst)])
                %data(i).data = d;
            end
        end
    end
    
end

