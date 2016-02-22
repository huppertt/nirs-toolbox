classdef FieldTrip
    
    properties
        mesh; 
        probe;
        prop;
    end
    
    methods
        %% Constructor
        function obj = FieldTrip( mesh, prop, probe)
            if nargin > 0, obj.mesh = mesh; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
        end
        
        %% Methods
        meas        = measurement( obj );
        [J,meas]    = jacobian( obj, type );
        mesh        = getFieldTripCFG( obj, type )
        
    end
    

        
end

