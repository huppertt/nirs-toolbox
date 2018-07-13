function list = pipelineToList( pipeline )
%% pipelineToList - Disassembles a pipeline into a cell array of modules.
%
% Args:
%     pipeline - an pipeline consisting of modules

    list = {};
    while ~isempty( pipeline )
        % get module
        m = pipeline;
        
        % strip off the pipeline
        m.prevJob = [];
        
        % append module
        list{end+1,1} = m;
        
        % get the rest of the pipeline
        pipeline = pipeline.prevJob;
    end
 
    % reverse list since we took 
    % modules from the end
    list = flipud(list);
    
end
