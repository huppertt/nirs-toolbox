classdef FunctionalData < nirs.Data
    %DATA Object to hold nirs data
    
    properties
        stimulus = {};        % struct containing stim vectors (vectors, names, types)
        demographics = table(); 	% table containing demographics (names, values)
    end

    methods
        %% Constructor
        function obj = FunctionalData( data, time, probe, Fm, stimulus, demographics )%, description )
            if nargin > 0, obj.data         = data;         end
            if nargin > 1, obj.time         = time;         end
            if nargin > 2, obj.probe        = probe;        end
            if nargin > 3, obj.Fm           = Fm;           end
            if nargin > 4, obj.stimulus     = stimulus;     end
            if nargin > 5, obj.demographics = demographics; end
        end
        
        %% Set/Get
        function obj = set.stimulus( obj, stim )
           assert( iscell(stim) )
           obj.stimulus = stim;
        end
        
        function obj = set.demographics( obj, demo )
           assert( istable( demo ) )
           obj.demographics = demo;
        end
    end
end

