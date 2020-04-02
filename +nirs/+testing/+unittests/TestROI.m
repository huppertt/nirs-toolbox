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
            obj.ROIs= table(NaN,NaN,{'ALL'},'VariableNames',{'source','detector','name'});

            % Test case #1
            % hbo- 1. low beta, high var (t=1), 2. high beta, low var (t=6)
            % hbr- 1. low beta, low var (t=-3),  2. high beta, high var (t=-2)
            channel_stats = nirs.core.ChannelStats();
            channel_stats.probe = probe;
            channel_stats.variables = [probe.link table(repmat({'Test1'},height(probe.link),1),'VariableNames',{'cond'})];
            channel_stats.dfe = 40;
            channel_stats.beta = [ 3 -3  6 -6 ]';
            channel_stats.covb = [ 9 -2  2 -3 ;
                                  -2  3 -1  2 ;
                                   2 -1  3 -2 ;
                                  -3  2 -2  9 ];
            
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
            hbo_inds = strcmp(obj.data1.variables.type,'hbo');
            hbr_inds = strcmp(obj.data1.variables.type,'hbr');
            true_betas(1,1) = sum(obj.data1.beta(hbo_inds)) / sum(hbo_inds);
            true_betas(2,1) = sum(obj.data1.beta(hbr_inds)) / sum(hbr_inds);
            true_covb(1,1) = mean(reshape(obj.data1.covb(hbo_inds,hbo_inds),[],1));
            true_covb(1,2) = mean(reshape(obj.data1.covb(hbo_inds,hbr_inds),[],1));
            true_covb(2,1) = mean(reshape(obj.data1.covb(hbr_inds,hbo_inds),[],1));
            true_covb(2,2) = mean(reshape(obj.data1.covb(hbr_inds,hbr_inds),[],1));
            
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
            weights = 1./diag(obj.data1.covb);
            weights2 = weights * weights';
            hbo_inds = strcmp(obj.data1.variables.type,'hbo');
            hbr_inds = strcmp(obj.data1.variables.type,'hbr');
            true_betas(1,1) = sum(weights(hbo_inds) .* obj.data1.beta(hbo_inds)) / sum(weights(hbo_inds));
            true_betas(2,1) = sum(weights(hbr_inds) .* obj.data1.beta(hbr_inds)) / sum(weights(hbr_inds));
            true_covb(1,1) = sum(sum(weights2(hbo_inds,hbo_inds) .* obj.data1.covb(hbo_inds,hbo_inds))) / sum(sum(weights2(hbo_inds,hbo_inds)));
            true_covb(1,2) = sum(sum(weights2(hbo_inds,hbr_inds) .* obj.data1.covb(hbo_inds,hbr_inds))) / sum(sum(weights2(hbo_inds,hbr_inds)));
            true_covb(2,1) = sum(sum(weights2(hbr_inds,hbo_inds) .* obj.data1.covb(hbr_inds,hbo_inds))) / sum(sum(weights2(hbr_inds,hbo_inds)));
            true_covb(2,2) = sum(sum(weights2(hbr_inds,hbr_inds) .* obj.data1.covb(hbr_inds,hbr_inds))) / sum(sum(weights2(hbr_inds,hbr_inds)));
                                    
            % Compare
            obj.verifyEqual(res.beta,true_betas);
            obj.verifyLessThan(sum(res.covb(:)-true_covb(:)),10^-12);
            
        end

    end
    
end

