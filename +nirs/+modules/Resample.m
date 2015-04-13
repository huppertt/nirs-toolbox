classdef Resample < nirs.functional.AbstractModule
  
    properties
        Fs = 4;
    end
    
    methods

        function obj = Resample( prevJob )
           obj.name = 'Resample';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                
                assert( obj.Fs < data(i).Fs )
                
                % resample data
                d = data(i).data;
                t = data(i).time;
                
                N = floor(t(end) * obj.Fs);
                new_t = t(1) + (0:N-1)' / obj.Fs;
                
                % anti-aliasing filter
                ord = floor( length(t) / 10 );
                Fc = obj.Fs/data(i).Fs;

                b = fir1(ord,Fc);
                
                for j = 1:2 % backward then forward
                    d = flipud( d );
                    d1 = d(1,:);
                    d = bsxfun(@minus,d,d1);
                    d = filter(b,1,d);
                    d = bsxfun(@plus,d,d1);
                end
                
                % interpolation
                d = interp1(t,d,new_t,'pchip');

                data(i).data = d;
                data(i).time = new_t;
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

