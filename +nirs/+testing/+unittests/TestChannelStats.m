classdef TestChannelStats < matlab.unittest.TestCase

    properties
        stats
    end
    
    methods
        function obj = TestChannelStats()
            % make a probe
            srcPos = [0 0 0];
            
            detPos = [-20 0 0;
                       20 0 0];

            link = [1 1 690;
                1 1 830;
                1 2 690;
                1 2 830];

            link = table(link(:,1), link(:,2), link(:,3), ...
                'VariableNames', {'source', 'detector', 'type'});

            probe = nirs.core.Probe(srcPos, detPos, link);
            
            % make channel stats object
            stats = nirs.core.ChannelStats();
            
            beta = [ones(4, 1); zeros(4,1)];
            covb = toeplitz([1 0.5 zeros(1, 6)]);
            
            dfe = 100;
            
            cond = {'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B'}';
            vars = [[link; link] table(cond)];
            
            stats.beta = beta;
            stats.covb = covb;
            stats.dfe = dfe;
            stats.variables = vars;
            stats.probe = probe;
            
            obj.stats = stats;
        end
    end
    
    methods (Test)
        function testConditions( obj )
            obj.verifyEqual(obj.stats.conditions, {'A', 'B'}');
        end
        
        function testTstats( obj )
            b = obj.stats.beta;
            covb = obj.stats.covb;
            t = b ./ sqrt(diag(covb));
            obj.verifyEqual(obj.stats.tstat, t);
        end

        function testPvals( obj )
            b = obj.stats.beta;
            covb = obj.stats.covb;
            t = b ./ sqrt(diag(covb));
            p = 2*tcdf(-abs(t), obj.stats.dfe);
        end

        function testQvals( obj )

        end
    end
end

