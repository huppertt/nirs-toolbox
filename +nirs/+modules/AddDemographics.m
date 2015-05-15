classdef AddDemographics < nirs.modules.AbstractModule
  
    properties
        demoTable
        varToMatch = 'subject';
    end
    
    methods

        function obj = AddDemographics( prevJob )
           obj.name = 'Add Demographics From Table';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            % columns of the table that arent varToMatch
            colNames = obj.demoTable.Properties.VariableNames;
            lst = find( ~strcmp( obj.varToMatch, colNames ) );
                
            for i = 1:length(data)

                % row idx of demo table
                idx = find( strcmpi( obj.demoTable.(obj.varToMatch), ...
                    data(i).demographics(obj.varToMatch) ) );
                
                
                % add demo info from row
                for j = 1:length(lst)
                    key = colNames{lst(j)};
                    val = obj.demoTable.(colNames{lst(j)}) (idx);
                    
                    if iscell( val )
                        val = val{1};
                    end
                    
                    data(i).demographics( key ) = val;
                end
                
            end
        end
    end
    
end

