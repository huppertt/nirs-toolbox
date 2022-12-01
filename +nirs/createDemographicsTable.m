function tbl = createDemographicsTable( data )
%% CREATEDEMOGRAPHCISTABLE - returns a table of demographics info per file
% 
% Args:
%     data - a list of nirs.core.Data or nirs.core.ChannelStats objects
%     
% Returns:
%     tbl  - a table containing demographics variables per item in data
    

keys=[]; 
for i=1:length(data)
   keys=[keys(:); data(i).demographics.keys(:)]; 
end;
keys=unique(keys);

    for j=1:length(keys)
        for i=1:length(data)
             if(ismember(keys(j),data(i).demographics.keys))
                 a=data(i).demographics(keys(j));
             end
        end
        if(isnumeric(a)); a=NaN; end;
        if(iscellstr(a)); a={''}; end;
        if(isstr(a)); a=''; end;
        if(ischar(a)); a=''; end;
     for i=1:length(data)
            
        if(~ismember(keys(j),data(i).demographics.keys))
            data(i).demographics(keys(j))=a;
        end
    end
end

    tbl=struct;    
    % loop over data files
    for i = 1:length(data)
       tbl(i)=tbl(1);
       flds=fields(tbl(i));
       for f=1:length(flds)
        tbl(i).(flds{f})=[];
       end
        % get demographics
        clear demo;
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
                        elseif(~isequaln(tbl(i).(demo{k}.keys{j}), demo{k}.values{j}))
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
                    keys=demo.keys{j};
                    keys(strfind(keys,'-'))='_';
                    keys(strfind(keys,' '))='_';
                    
                    if(isa(demo.values{j},'Dictionary'))
                        tbl(i).(keys) = {demo.values{j}};
                    elseif(isa(demo.values{j},'char'))
                        tbl(i).(keys) = cellstr(demo.values{j});
                    else
                    tbl(i).(keys) = demo.values{j};
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