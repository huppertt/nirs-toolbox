classdef KalmanARWLS < handle
    properties( SetAccess = private )
        modelKF;
        arKF;
    end
    
    properties( Dependent )
        B;
        covB;
        a;
        tstat;
        dfe;
    end
    
    properties( Access = private )
        ybuffer;
        Xbuffer;
    end
    
    methods
        % boa constructer
        function obj = KalmanARWLS( modelKF, arKF )
            obj.modelKF = modelKF;
            obj.arKF = arKF;
            obj.ybuffer = zeros(obj.arKF.Pmax,1);
            obj.Xbuffer = zeros(obj.arKF.Pmax,obj.modelKF.np);
        end
        
        % getters
        function B = get.B(obj)
            B = obj.modelKF.B; 
        end
        function a = get.a(obj)
            a = obj.arKF.a; 
        end
        function covB = get.covB(obj)
            covB = obj.modelKF.P;
        end
        function tstat = get.tstat(obj)
            tstat = obj.B ./ sqrt(diag( obj.covB ));
        end
        function dfe = get.dfe(obj)
            dfe = obj.modelKF.n - obj.modelKF.np;
        end
        
        % update method
        function update(obj,y,X,F,u)
            if nargin < 4
                F = 1;
            end
            if nargin < 5
                u = 0;
            end
            
            % whiten with ar coefs
            a = obj.arKF.a(2:end);
            f = [1; -a];
            
            yf = f' * [y; obj.ybuffer(1:length(a))];
            Xf = sum([X; obj.Xbuffer(1:length(a),:)] .* repmat(f,[1 obj.modelKF.np]));
            
            % update model
            obj.modelKF.update(yf,Xf,F,u);
            
%             % update ar
%          	obj.arKF.update(...
%                 y-X*obj.modelKF.B ...
%                 );
            
            % update ar
         	obj.arKF.update(...
                [y; obj.ybuffer]-[X; obj.Xbuffer]*obj.modelKF.B ...
                );
            
            obj.updateBuffers(y,X);
        end
    end
    
    methods( Access = private )
        function updateBuffers(obj,y,X)
            obj.ybuffer = [y; obj.ybuffer(1:end-1)];
            obj.Xbuffer = [X; obj.Xbuffer(1:end-1,:)];
        end
        
    end
    
end

