classdef WaveletImageRecon
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X
        W
        C
    end
    
    properties( SetAccess = private )
        U
        S
        V
        L
    end
    
    properties( Access = private )
        flag = true;
    end
    
    methods
        function obj = WaveletImageRecon( X, W )
            if nargin > 0, obj.X = X; end
            if nargin > 1, obj.W = W; end
        end
        
        function obj = set.X( obj, X )
           obj.X = X;
           obj.flag = true;
        end
        
        function obj = set.W( obj, W )
           obj.W = W;
           obj.flag = true;
        end
        
        function obj = set.X( obj, X )
           obj.X = X;
           obj.flag = true;
        end
        
        function b = solve( obj, y )
            if obj.flag == true
               obj = obj.calculateIntermediates(); 
            end
            
            iW = obj.V*((obj.W*obj.V)\eye(size(obj.V,2)));
            XW = obj.X*iW;
            
            iXW = lambda*XW'*pinv(XW*XW'*lambda + eye(size(XW,1)));
            what = iXW*y;

            H = XW*iXW;

            dfe = trace(H'*H)^2 / trace(H'*H*H*H');

            lambda = (what'*what)/dfe;
        end
        
        function obj = obj.calculateIntermediates( obj )
            if isempty 
        end
    end
    
end