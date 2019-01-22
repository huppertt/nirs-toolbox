function tbl = createDemographicsTable( data )
%% CREATEDEMOGRAPHCISTABLE - returns a table of demographics info per file
% 
% Args:
%     data - a list of nirs.core.Data or nirs.core.ChannelStats objects
%     
% Returns:
%     tbl  - a table containing demographics variables per item in data
    
    tbl=struct;    
    % loop over data files
    for i = 1:length(data)
       tbl(i)=tbl(1);
       flds=fields(tbl(i));
       for f=1:length(flds)
        tbl(i).(flds{f})=[];
       end
        % get demographics
        demo = data(i).demographics;
        if(iscell(demo))
           
            % case of hyperscan data
            for k=1:length(demo)
               
                if(istable(demo{k}))
                    tbl(i)=demo{k};
                   
                end
                % create struct
                if(isprop(demo{k},'keys'))
                    for j = 1:length(demo{k}.keys)
                        if(~isfield(tbl(i),demo{k}.keys{j}) | isempty(demo{k}.values{j}))
                            tbl(i).(demo{k}.keys{j}) = demo{k}.values{j};
                        elseif(isempty(tbl(i).(demo{k}.keys{j})))
                            tbl(i).(demo{k}.keys{j}) = demo{k}.values{j};
                        elseif(isstr(demo{k}.values{j}) && ~strcmp(tbl(i).(demo{k}.keys{j}),demo{k}.values{j}))
                            
                            tbl(i).([demo{k}.keys{j} '_' num2str(k)]) = demo{k}.values{j};
                        elseif(tbl(i).(demo{k}.keys{j}) ~= demo{k}.values{j})
                            tbl(i).([demo{k}.keys{j} '_' num2str(k)]) = demo{k}.values{j};
                        end
                    end
                end
            end
        else
            if(istable(demo))
                tbl=demo;
                return
            end
            % create struct
            if(isprop(demo,'keys'))
                for j = 1:length(demo.keys)
                    if(isa(demo.values{j},'Dictionary'))
                         tbl(i).(demo.keys{j}) = {demo.values{j}};
                    elseif(isa(demo.values{j},'char'))
                        tbl(i).(demo.keys{j}) = cellstr(demo.values{j});
                    else
                    tbl(i).(demo.keys{j}) = demo.values{j};
                    end
                end
            end
        end
    end
    % covert to table
    if(exist('tbl'))
        tbl = struct2table(tbl);
    else
        tbl=table;
    end
end