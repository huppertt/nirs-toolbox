function pipeline = listToPipeline( list )
%% listToPipeline - Assembles a cell array of modules into a single pipline.
%
% Args:
%     list - cell array of modules
    
    pipeline = [];
    for i = 1:length(list)
        next = list{i};
        next.prevJob = pipeline;
        
        pipeline = next;
    end
    
end

