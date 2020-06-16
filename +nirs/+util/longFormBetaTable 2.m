function tbl = longFormBetaTable( S )
%% longFormBetaTable - returns a table with all regression stats
% 
% Args:
%     S - list of ChannelStats or ROI stats from first-level analysis
    if isa(S, 'nirs.core.ChannelStats')
        tbl = tableFromChannelStats( S )
    else
        tbl = tableFromTables(S);
    end
    
end

function tbl = tableFromChannelStats( S )
    t = S(1).table();
    
    for i = 1:size(t,1)
        colNames{i} = strjoin({ ...
            ['S' num2str(t.source(i))], ...
            ['D' num2str(t.detector(i))], ...
            t.type{i}, ...
            t.cond{i}}, '_');
    end
    
    for i = 1:length(S)
        beta(i,:) = S(i).beta';
    end

    tbl = array2table(beta);
    tbl.Properties.VariableNames = colNames;
    
    tbl = [nirs.createDemographicsTable(S) tbl];
end

function tbl = tableFromTables( S )
    beta = [];
    for i = 1:length(S)
        beta(i,:) = S{i}.Beta;
    end

    for i = 1:size(S{1},1)
        colNames{i} = strjoin( {S{1}.ROI{i} S{1}.Contrast{i}}, '_' );
    end

    tbl = array2table(beta);
    tbl.Properties.VariableNames = colNames;

end