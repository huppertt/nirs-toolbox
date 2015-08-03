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
            
            beta = [1; 2; 3; 4; zeros(4,1)];
            covb = toeplitz([1 0.5 zeros(1, 6)]);
            
            dfe = 100;
            
            cond = {'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B'}';
            vars = [[link; link] table(cond)];
            
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type'});
            
            stats.beta  = beta(idx);
            stats.covb  = covb(idx, idx);
            stats.dfe   = dfe;
            
            stats.variables = vars;
            stats.probe    	= probe;
            
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
            obj.verifyEqual(obj.stats.p, p);
        end

        function testQvals( obj )
            % it is possible that q can be greater than p
            % but if this test fails, verify fdr code
            obj.verifyLessThanOrEqual( obj.stats.p, obj.stats.q )
        end
        
        function testGetCritT( obj )
            tcrit = -tinv(0.025, obj.stats.dfe);
            obj.verifyEqual( obj.stats.getCritT('p < 0.05'), tcrit)
        end
        
        % most important tests
        function testTtest( obj )
            c = [1 1; 1 -1]; % A+B; A-B
            s = obj.stats.ttest(c);
            
            % manual calculation
            C = kron(eye(4), c);
            b = C*obj.stats.beta;
            covb = C*obj.stats.covb*C';
            
            obj.verifyEqual(s.beta, b)
            obj.verifyEqual(s.covb, covb)
            
            % make sure ordering of variables wont effect results
            newObj = obj.stats;
            idx = [1 3 2 8 7 4 5 6]';
            
            newObj.beta = obj.stats.beta(idx);
            newObj.covb = obj.stats.covb(idx, idx);
            newObj.variables = obj.stats.variables(idx, :);
            
            newS = newObj.ttest(c);
            
            obj.verifyEqual(s, newS);
        end
        
        function testFtest( obj )
            m = eye(2) > 0;
            
            % check simple relation of F to t^2
            f = obj.stats.ftest(m);
            
            obj.verifyEqual(f.F, obj.stats.tstat.^2);
            obj.verifyEqual(f.p, obj.stats.p);
            
            % more complicated
            m = [1 0; 1 1] > 0;
            f = obj.stats.ftest(m);
            
            M = kron(eye(4), m) > 0;
            
            b = obj.stats.beta;
            covb = obj.stats.covb;
            
            for i = 1:size(M,1)
                idx = M(i,:);
                T2(i,1)     = b(idx)'*inv(covb(idx,idx))*b(idx);
                df1(i,1)    = sum(idx);
                df2(i,1)    = obj.stats.dfe - df1(i) + 1;
                F(i,1)      = df1(i) / df2(i) / obj.stats.dfe * T2(i);
                p(i,1)      = fcdf(1./F(i), df2(i), df1(i));
            end
            
            
            
        end
        
        
        function testJointTest( obj )
            
        end
    end
end

