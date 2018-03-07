classdef AuxToStim < nirs.modules.AbstractModule
%%  Create stimulus events from auxillary data

    properties
        replace_existing = false;
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
                
                [stim_vectors, stim_names] = nirs.util.aux2stim_kmeans( aux.data );
                
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
