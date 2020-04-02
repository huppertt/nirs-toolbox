classdef HRFwienerfilter < nirs.modules.AbstractModule
%% Implementaton of a wiener filter using the HRF kernel
% This is not exactly a wiener filter, TODO- fix it
    properties
     
    end

    methods

        function obj =  HRFwienerfilter( prevJob )
           obj.name = 'HRF Convolve';
         if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                Fs = round(data(i).Fs);
                d = data(i).data;
                hrf=convert(nirs.design.basis.Canonical,[1; zeros(20*Fs,1)],[0 1/Fs]);
                f=[1; -hrf];
                
                d=filtfilt(hrf,1,filtfilt(f,1,d));
                data(i).data = d;
                
                               
            end
        end
    end
    
end

