classdef HRFconvolve < nirs.modules.AbstractModule
%% COnvolves an EEG signal with the HRF kernel
%
% Options:
%     ncomp - % number of components to remove
    properties
     
    end

    methods

        function obj = BandPassFilter( prevJob )
           obj.name = 'HRF Convolve';
         if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                 hrf=convert(nirs.design.basis.Canonical,[1; zeros(20*data(i).Fs,1)],[0 1/data(i).Fs]);
                d=filter(hrf,1,d);
                
                data(i).data = d; 
                
                               
            end
        end
    end
    
end

