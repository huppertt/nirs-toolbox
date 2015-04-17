classdef Data
    %DATA Object to hold nirs data
    
    properties
        description;
        data;               % channel time series in columns
        probe;              % object describing geometry
        time;               % vector of time points
        Fm = 0;             % modulation frequency in MHz: 0 for CW; 110 for ISS
        
    end
    
    properties( Dependent = true )
        Fs = 0;             % sampling frequency in Hz
    end
    
    properties( Access = private )
        
    end
    
    methods
        %% Constructor
        function obj = Data( data, time, probe, Fm )%, description )
            if nargin > 0, obj.data = data; end
            if nargin > 1, obj.time = time; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Set/Get
        function obj = set.Fm( obj, Fm )
            assert( isscalar(Fm) && obj.Fm >= 0 )
            obj.Fm = Fm;
        end

        function obj = set.data( obj, data )
            assert( ismatrix(data) )
            obj.data = data;
        end
        
        function obj = set.probe( obj, probe )
            assert( isobject(probe) || isstruct(probe) )
            obj.probe = probe;
        end
        
        function obj = set.time( obj, time )
           assert( isvector(time) )
           obj.time = time(:);
        end
        
        function out = get.Fs( obj )
            if length(obj.time) > 1
                out = 1 / ( obj.time(2) - obj.time(1) );
            else
                out = NaN;
            end
        end
        
        % show data
        function draw( obj, lstChannels )
            if nargin == 1
                plot(obj.data)
            else
                plot( obj.data(:,lstChannels) ) 
            end
        end
        
    end
end

