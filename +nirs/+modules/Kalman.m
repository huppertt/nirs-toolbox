classdef Kalman < nirs.modules.AbstractModule
    %% RTS-Kalman filter of data
    %
    % Options:
    %     Q - process noise
    
    
    properties
        Q = 0; % process  noise
        verbose=true;
    end
  
    methods
        function obj = Kalman( prevJob )
            obj.name = 'Kalman';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            
           
            
            for i = 1:numel(data)
                if(obj.verbose)
                    disp(['running ' data(i).description ' (' num2str(i) ' of ' num2str(numel(data)) ')']);
                end
                if(~isempty(data(i).data))
                    try
                        data(i).data=nirs.math.kalman_rts(data(i).data',[],obj.Q)';
                    catch
                        warning(['problem with file: ' num2str(i) ' ' data(i).description]);
                    end
                end
            end
        end
    end
    
end

