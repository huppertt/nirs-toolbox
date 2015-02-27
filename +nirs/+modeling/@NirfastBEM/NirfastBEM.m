classdef NirfastBEM
    
    properties
        mesh; 
        probe;
        prop;

        Fm = 0;
    end
    
%     properties(SetAccess = private)
%         nLayers;
%     end
    
    methods
        %% Constructor
        function obj = NirfastBEM( mesh, prop, probe, Fm )
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
        mesh = getNirfastMeshes( obj )
        [J,meas] = jacobian( obj );
        [J,meas] = layeredJacobian( obj );
        
        
%         [J,meas] = spectralJacobian( obj );
        
    end
    

        
end

