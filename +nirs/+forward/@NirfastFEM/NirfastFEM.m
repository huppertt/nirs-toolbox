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
        
        %% Methods
        meas        = measurement( obj );
        [J,meas]    = jacobian( obj, type );
        mesh        = getNirfastMeshes( obj, type )
        
    end
    

        
end

