classdef AddDemographics < nirs.modules.AbstractModule
%% AddDemographics - Adds demographics info to data given a table.
% 
% Options: 
%     demoTable     - a matlab table with the demographics info
%     varToMatch    - the column of demoTable which also matches an
%                     existing demographics variable (e.g. match on subject
%                     name)
%     allowMissing  - Flag to prevent erroring out on missing subjects (default: false)
    properties
        demoTable                  % a matlab table with the demographics info
        varToMatch = 'subject';    % the column of demoTable to match files with
        allowMissing = false;
    end
    
    methods

        function obj = AddDemographics( prevJob )
           obj.name = 'Add Demographics From Table';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
              %remove a few bad rows (bad excel format)
            f=obj.demoTable.Properties.VariableNames;
            for i=1:length(f)
                if(~isempty(strfind(f{i},'Var')))
                    obj.demoTable.(f{i})=[];
                end
            end

            
            % columns of the table that arent varToMatch
            colNames = obj.demoTable.Properties.VariableNames;
            if ~iscell(obj.varToMatch)
                lst = find( ~strcmp( obj.varToMatch, colNames ) );
            else
                lst = find(~ismember(colNames,obj.varToMatch));
            end

            idx=[];
            varToMatchA_idx=find(strcmpi(obj.demoTable.Properties.VariableNames,obj.varToMatch));
            if(isempty(varToMatchA_idx))
                if obj.allowMissing
                    warning('All entries missing in demographics table for %s!',obj.varToMatch);
                else
                    error('All entries missing in demographics table for %s!',obj.varToMatch);
                end
            end
            
            
            for i = 1:numel(data)
                    

                if ~iscell(obj.varToMatch)
                    % Make this case-insensitive
                    if(~isempty(varToMatchA_idx))
                        varToMatchA=obj.demoTable.Properties.VariableNames{varToMatchA_idx};
                         % row idx of demo table
    
                         if(isstr(data(i).demographics(obj.varToMatch)) | iscellstr(data(i).demographics(obj.varToMatch)))
                            idx = find( strcmpi( obj.demoTable.(varToMatchA), data(i).demographics(obj.varToMatch) ) );
                         else
                             idx=find(obj.demoTable.(varToMatchA)==data(i).demographics(obj.varToMatch) );
                         end
                    end
                     
                else
                    if(~isempty(varToMatchA_idx))
                        varToMatchA=obj.demoTable.Properties.VariableNames(find(ismember(lower(obj.demoTable.Properties.VariableNames),lower(obj.varToMatch))));
                        % row idx of demo table
                        T = cell2table(data(i).demographics(obj.varToMatch),'VariableNames',varToMatchA);
                      %  [~,idx] = ismember(T,obj.demoTable(:,varToMatchA));
                        idx = find(ismember(obj.demoTable(:,varToMatchA),T));
                    end
                end
                
               
                
                if(isempty(idx)&&~isempty(varToMatchA_idx))
                    if obj.allowMissing
                        warning(['Missing entry: ' data(i).demographics(obj.varToMatch)]);
                    else
                        error(['Missing entry: ' data(i).demographics(obj.varToMatch)]);
                    end
                end
                if(numel(idx)>1)
                    warning(['Duplicate entry: ' data(i).demographics(obj.varToMatch)]);
                end
                    
                % add demo info from row
                for j = 1:length(lst)
                    key = colNames{lst(j)};
                    coldata = obj.demoTable.(colNames{lst(j)});
                    
                    if isempty(idx)
                        if isnumeric(coldata)
                            val = nan;
                        else
                            val = '';
                        end
                    else
                        val = coldata(idx);

                        if iscell( val )
                            val = val{1};
                        end
                    end
                    
                    data(i).demographics( key ) = val;
                end
                
            end
        end
    end
    
end

