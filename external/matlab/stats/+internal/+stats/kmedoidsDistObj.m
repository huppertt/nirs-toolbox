classdef kmedoidsDistObj
%KMEDOIDSDISTOBJ - an internal helper object to abstract distance 
%calculations in kmedoids
%
% Usage
%
% distObj = kmedoidsDistObj(X,distance);
%
% where distance is either 'sqEuclidean' or a distance metric supported by
% pdist/pdist2.
%
% out = distObj.pdist(Y);
%
% performs pdist on Y using the distance metric
%
% kmedoidsDistObj is not intended to be called by users and may change or
% be removed in a future release.

% Copyright MathWorks 2014

    properties (SetAccess = private)
        distance
        standardEucScaling
        mahalScaling
    end
    
    methods
        
        function distObj = kmedoidsDistObj(X,distance)
            distObj.distance = distance;
            % better validation needed, or completed in kmedoids?
            if ischar(distObj.distance)
                switch distObj.distance
                    case 'seuclidean'
                        distObj.standardEucScaling = nanstd(X);
                        distObj.mahalScaling       = [];
                    case 'mahalanobis'
                        distObj.standardEucScaling = [];
                        distObj.mahalScaling       = nancov(X);
                    otherwise
                        distObj.standardEucScaling = [];
                        distObj.mahalScaling       = [];
                end
            else
                distObj.standardEucScaling = [];
                distObj.mahalScaling       = [];
            end
        end
        
        function out = pdist(distObj,X)
            if ischar(distObj.distance)
                switch distObj.distance
                    case 'seuclidean'
                        out = pdist(X,'seuclidean',distObj.standardEucScaling);
                    case 'mahalanobis'
                        out = pdist(X,'mahalanobis',distObj.mahalScaling);
                    case 'sqeuclidean'
                        out = pdist(X,'euclidean');
                        out = out.^2;
                    otherwise
                        out = pdist(X,distObj.distance);
                end
            else
                % function handle
                out = pdist(X,distObj.distance);
                
            end
        end
        
        function out = pdist2(distObj,X1,X2)
            if ischar(distObj.distance)
                switch distObj.distance
                    case 'seuclidean'
                        out = pdist2(X1,X2,'seuclidean',distObj.standardEucScaling);
                    case 'mahalanobis'
                        out = pdist2(X1,X2,'mahalanobis',distObj.mahalScaling);
                    case 'sqeuclidean'
                        out = pdist2(X1,X2,'euclidean');
                        out = out.^2;
                    otherwise
                        out = pdist2(X1,X2,distObj.distance);
                end
            else
                out = pdist2(X1,X2,distObj.distance);
            end
        end
        
    end
    
end