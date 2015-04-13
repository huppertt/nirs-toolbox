classdef FunData < nirs.Data
    %DATA Object to hold nirs data
    
    properties
        stimulus        = HashTable();	% struct containing stim vectors (vectors, names, types)
        demographics    = HashTable();	% table containing demographics (names, values)
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
           assert( isa(stim,'HashTable') )
           obj.stimulus = stim;
        end
        
        function obj = set.demographics( obj, demo )
           assert( isa(demo,'HashTable') )
           obj.demographics = demo;
        end
        
        %% show data
        function draw( obj, lstChannels )
            % get data
            if nargin == 1
                lstChannels = 1:size(obj.data,2);
            end
            
            t = obj.time;
            d = obj.data(:,lstChannels);
            
            % get stim vecs
            s = []; k = obj.stimulus.keys;
            for i = 1:length( obj.stimulus.keys )
                s = [s obj.stimulus.values{i}.getStimVector( t )];
            end

            % plots
            gca, hold on
            
            % data min/max/size
            dmax = max( d(:) );
            dmin = min( d(:) );
            dsize = (dmax-dmin);
            
            % plot stim blocks if available
            if ~isempty(s) 
                % min/max of axes
                pmin = dmin - 0.3*dsize;
                pmax = dmax + 0.2*dsize;

                % adjust amplitude so stims are visible
                s = 0.15*dsize*s + dmin - 0.25*dsize;
                
                % plot
                plot( t, s, 'LineWidth', 3 )
                
                % legend
                l = legend(k{:});
                set(l,'Interpreter', 'none')
            else
                % min/max of axes
               	pmin = dmin - 0.1*dsize;
                pmax = dmax + 0.1*dsize;
            end
            
            % plot data
            plot( t, d ) 
            
            % axes limits
            xlim( [min(t) max(t)] )
            ylim( [pmin pmax] )

        end
        
    end
end

