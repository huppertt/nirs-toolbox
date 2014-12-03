classdef Data
    %DATA Object to hold nirs data
    
    properties
        data;               % channel time series in columns
        probe;              % object describing geometry
        Fs = 0;           	% sampling frequency in Hz
        Fm = 0;             % modulation frequency in MHz: 0 for CW; 110 for ISS
        description;        
    end
    
    methods
        %% Constructor
        function obj = Data( data, probe, Fs, Fm, description )
            if nargin > 0, obj.data = data; end
            if nargin > 1, obj.probe = probe; end
            if nargin > 2, obj.Fs = Fs; end
            if nargin > 3, obj.Fm = Fm; end
            if nargin > 4, obj.description = description; end
        end
        
        %% Set/Get
        function obj = set.Fs( obj, Fs )
            assert( isscalar(Fs) && obj.Fs >= 0 )
            obj.Fs = Fs;
        end
        
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
        
        function obj = set.description( obj, description )
            assert( ischar( description ) || isempty( description ) )
          	obj.description = description;
        end
        
    end
end

