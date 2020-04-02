classdef RemoveShortScans< nirs.modules.AbstractModule
%% RemoveShortScan - Removes files shorter then minimum time in length
% 
    properties
        mintime=30;  % time in seconds
        
    end
    methods
        function obj = RemoveShortScans( prevJob )
           obj.name = 'Remove Short Scans';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            lstbad=[];
            for i = 1:numel(data)
              if((max(data(i).time)-min(data(i).time))<obj.mintime)
                  lstbad=[lstbad; i];
              end
            end
            data(lstbad)=[];
        end
    end
    
end

