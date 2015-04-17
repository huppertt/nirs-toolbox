function pipeline = listToPipeline( list )

    pipeline = [];
    for i = 1:length(list)
        next = list{i};
        next.prevJob = pipeline;
        
        pipeline = next;
    end
    
end

