classdef ForwardModel
    %FWDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    methods( Abstract )
        [J,meas] = jacobian( obj );
        meas = measurement( obj );
    end
    
end

