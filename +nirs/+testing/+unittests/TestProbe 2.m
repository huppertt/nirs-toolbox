classdef TestProbe < matlab.unittest.TestCase

    properties
        probe
    end
    
    methods
        function obj = TestProbe()
            srcPos = [0 0 0];
            
            detPos = [-20 0 0;
                       20 0 0];

            link = [1 1 690;
                1 1 830;
                1 2 690;
                1 2 830];

            link = table(link(:,1), link(:,2), link(:,3), ...
                'VariableNames', {'source', 'detector', 'type'});

            obj.probe = nirs.core.Probe(srcPos, detPos, link);
        end
    end
    
    methods (Test)
        function testDistances( obj )
            obj.verifyEqual(obj.probe.distances, 20*ones(4,1));
        end
        
        function testSwapSD( obj )
           swapped = obj.probe.swapSD();
           obj.verifyEqual(obj.probe.srcPos, swapped.detPos);
           obj.verifyEqual(obj.probe.detPos, swapped.srcPos);
           obj.verifyEqual(obj.probe.distances, swapped.distances);
        end
    end
    
end

