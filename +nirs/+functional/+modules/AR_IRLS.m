classdef AR_IRLS < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        hpf_Fc = 1/125;
        basis = {};
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
                t = data(i).time;
                Fs = data(i).Fs;
                
                % generate design matrix
                stim = data(i).stimulus;
                
                X = []; names = {}
                for j = 1:length(stim)
                    
                    basis = obj.basis{j};
                    x = basis.convert( stim{j}.getStimVector( t ) );
                    
                    if size(x,2) > 1
                        for k = 1:size(x,2)
                            names{end+1} = [stim{j}.name '_' sprintf('%02i',k)];
                        end
                    else
                        names{end+1} = stim{j}.name;
                    end
                    
                    X = [X x];
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

