classdef AuxToStim < nirs.modules.AbstractModule
%%  Create stimulus events from auxillary data

    properties
        method = 'threshold'; % or 'kmeans'
        replace_existing = false;
        threshold = [];
        verbose = false;
    end
    
    methods

        function obj = AuxToStim( prevJob )
           obj.name = 'Create stimulus events from auxillary data';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            

            for i = 1:length(data)
                
                aux = data(i).auxillary('aux');
                if isempty(aux)
                    continue;
                end
                
                if obj.replace_existing
                    data(i).stimulus = Dictionary();
                end
                
                if obj.verbose
                    fprintf('Converting AUX to stim (%i/%i)\n',i,length(data));
                end
                
                if strcmpi(obj.method,'threshold')
                    [stim_vectors, stim_names] = nirs.util.aux2stim( aux.data , obj.threshold );
                elseif any(strcmpi({'kmeans','k-means'},obj.method))
                    [stim_vectors, stim_names] = nirs.util.aux2stim_kmeans( aux.data , obj.threshold , obj.verbose );
                else
                    error('Unknown method: %s',obj.method);
                end
                
                for j = 1:length(stim_names)
                    
                    stim_event = nirs.design.vector2event( aux.time , stim_vectors(:,j) , stim_names{j} );

                    if ~isempty(stim_event.onset)
                        data(i).stimulus(stim_names{j}) = stim_event;
                    end
                    
                end
                
            end
        end
    end
end
