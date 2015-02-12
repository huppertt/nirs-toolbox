classdef AR_IRLS < nirs.jobs.AbstractJob
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        hpfCutoff = 1/125;
        stimTypes = containers.Map;
        constant = true;
    end
    
    methods

        function obj = AR_IRLS( prevJob )
           obj.name = 'GLM via AR(P)-IRLS';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                Fs = data(i).Fs;
                
                % generate design matrix
                stim = data(i).stimulus;
                keys = stim.keys;
                
                X = []; names = {}
                for j = 1:length(keys)
                    names = {}
                    
                    converter = stimTypes( keys(j) );
                    
                    
                end
                
                
                % call ar_irls
                
                
                % output stats
                
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

