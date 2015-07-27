classdef NirfastBEM
    
    properties
        mesh; 
        probe;
        prop;

        Fm = 0;
    end
    
    methods
        %% Constructor
        function obj = NirfastBEM( mesh, prop, probe, Fm )
            if nargin > 0, obj.mesh = mesh; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Methods
        meas = measurement( obj );
        mesh = getNirfastMeshes( obj )
        [J,meas] = jacobian( obj, type );
        [J,meas] = layeredJacobian( obj );
        
    end
    

        
end

