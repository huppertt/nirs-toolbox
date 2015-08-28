classdef MCXLab
    %MCXFWDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties%( Constant )
        gpuId = 1;
    end
    
    properties
       probe;
       prop;
       image; 
        directory = [getenv('TMPDIR') filesep 'tmp'...
            filesep num2str(randi(2^32-1))];

        Fm = 100;
        
        nPhotons = 1e7;
        nTimeGates = 40;
        timeStep = 1/100e6/300;
        nRepetitions = 1;
    end
    
    properties(SetAccess = private)
        nLayers;
        cleanup = true;
    end
    
    methods
        %% Constructor
        function obj = MCXLab(image, prop, probe, Fm )
            if nargin > 0, obj.image = image; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Set/Get
        
        function obj = set.directory( obj, d )
            assert( ischar( d ) )
            obj.directory = d;
            obj.cleanup = false;
        end
        
        %% Methods
        saveFluence( obj );
        meas = measurement( obj );
        meas = timeResolvedMeas( obj );
        [J,meas] = jacobian( obj );
        cfg = getConfig( obj, idx, idxType );
        
    end
    

        
end

