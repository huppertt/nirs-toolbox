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
        
        function obj = set.probe(obj,probe)
            if(~isa(probe,'nirs.core.Probe1020'))
                warning('probe must be a 3D registered probe');
            end
            if(all(probe.optodes.Z==0) & ~isa(probe,'nirs.core.Probe'))
                disp('warning: changing probe to 3D using "swap_reg" function');
                probe=probe.swap_reg;
            end
            obj.probe=probe;
        end
        
        
    end
    

        
end

