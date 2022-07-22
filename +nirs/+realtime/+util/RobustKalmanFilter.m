classdef RobustKalmanFilter < nirs.realtime.util.KalmanFilter
%     properties
%         Q = 0; 	% process noise
%     end
%     
%     properties( SetAccess = private )
%         P;      % covariance matrix
%         B;      % parameters
%         E = 0; 	% measurement noise
%         n = 0;	% number of updates
%         np;     % number of parameters
%     end
    properties
        w;
        tune = 4.685;
    end
    
    methods
        function obj = RobustKalmanFilter( Q, tune )
            obj = obj@nirs.realtime.util.KalmanFilter( Q );
            if nargin > 1
                obj.tune = tune;
            end
        end
        
        function update(obj,y,X,F,u)
            if nargin < 4
                F = 1;
            end
            if nargin < 5
                u = 0;
            end
% for i = 1:length(ty)
%     y = ty(i);
%     X = tX(i,:);
            obj.n = obj.n + 1;
            
        	% update measurement noise
            ehat = y - X*obj.B;
            obj.updateE( ehat );
            
            % prediction step
            obj.prediction_step( F,u );

            % weights
            if obj.n > 40
                w = obj.bisquare( ehat,obj.E,obj.tune );
            else
                w = 1;
            end
            
            % update B and P
            obj.update_step( w*ehat,w*X );
            
            obj.w = w;
% end
        end
    end
    
    methods( Static )
        function w = bisquare( x,sig,c )
            w = ( (1-(x/sig/c).^2) ) .* ( abs(x/sig) < c );
        end
    end
end