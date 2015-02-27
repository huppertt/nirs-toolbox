classdef NirfastFEM
    
    properties
        mesh; 
        probe;
        prop;

        Fm = 0;
    end
    
    methods
        %% Constructor
        function obj = NirfastFEM( mesh, prop, probe, Fm )
            if nargin > 0, obj.mesh = mesh; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Set/Get
%         function obj = set.mesh( obj, image )
%             assert( isa( image.vol,'uint8' ) )
%             obj.image = image;
%             obj.nLayers = max( obj.image.vol(:) );
%         end
        
%         function obj = set.probe( obj, probe )
%             assert( probe.makeUniqueProbe().isValid() )
%             obj.probe = probe;
%         end
        
%         function obj = set.directory( obj, d )
%             assert( ischar( d ) )
%             obj.directory = d;
%             obj.cleanup = false;
%         end
        
        %% Methods
        meas = measurement( obj );
        [J,meas] = jacobian( obj, type );
        mesh = getNirfastMeshes( obj )
        
    end
    

        
end

