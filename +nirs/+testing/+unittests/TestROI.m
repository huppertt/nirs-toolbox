classdef TestROI < matlab.unittest.TestCase

    properties
        ROIs
        data1
    end
    
    methods
        function obj = TestROI()
            
            % Create probe
            srcPos = [0 0 0];
            
            detPos = [ 0 30 0;
                30 0 0];
            
            s = [1 1 1 1]';
            d = [1 1 2 2]';
            t = {'hbo','hbr','hbo','hbr'}';
            link = table(s,d,t,'VariableNames', {'source', 'detector', 'type'});
            probe = nirs.core.Probe(srcPos, detPos, link);
            
            % Create ROI
            obj.ROIs = table({},{},{},'VariableNames',{'source','detector','name'});
            obj.ROIs(1,:) = array2table({[1 1],[1 2],'ALL'});

            % Test case #1
            % hbo- 1. low beta, high var (t=1), 2. high beta, low var (t=6)
            % hbr- 1. low beta, low var (t=-3),  2. high beta, high var (t=-2)
            channel_stats = nirs.core.ChannelStats();
            channel_stats.probe = probe;
            channel_stats.variables = [probe.link table(repmat({'Test1'},height(probe.link),1),'VariableNames',{'cond'})];
            channel_stats.dfe = 40;
            channel_stats.beta = [ 3 -3  6 -6 ]';
            channel_stats.covb = [ 9  0  0  0 ;
                                   0  1  0  0 ;
                                   0  0  1  0 ;
                                   0  0  0  9 ];
            
            obj.data1 = channel_stats;
            
        end
    end
    
    methods (Test)
        
        function testROI1_unweighted( obj )

            % Module calculation
            job = nirs.modules.ApplyROI();
            job.listOfROIs = obj.ROIs;
            job.weighted = false;          
            res = job.run( obj.data1 );
            
            % Manual calculation
            true_betas = (obj.data1.beta([1 2]) + obj.data1.beta([3 4])) ./ 2;
            true_covb = (obj.data1.covb([1 2],[1 2]) + obj.data1.covb([3 4],[3 4])) ./ 2;
            
            % Compare
            obj.verifyEqual(res.beta,true_betas);
            obj.verifyEqual(res.covb,true_covb);
            
        end
        
        function testROI1_weighted( obj )

            % Module calculation
            job = nirs.modules.ApplyROI();
            job.listOfROIs = obj.ROIs;
            job.weighted = true;          
            res = job.run( obj.data1 );
            
            % Manual calculation
            weights = 1./sqrt(diag(obj.data1.covb));
            true_betas = (weights([1 2]) .* obj.data1.beta([1 2])...
                        + weights([3 4]) .* obj.data1.beta([3 4])) ./ (weights([1 2])+weights([3 4]));
            true_covb = diag((weights([1 2]) .* diag(obj.data1.covb([1 2],[1 2])) ... 
                            + weights([3 4]) .* diag(obj.data1.covb([3 4],[3 4]))) ./ (weights([1 2])+weights([3 4])));
            
            % Compare
            obj.verifyEqual(res.beta,true_betas);
            obj.verifyLessThan(sum(res.covb(:)-true_covb(:)),10^-12);
            
        end

    end
    
end

