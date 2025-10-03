function queue = createParallelProgressBar(totalIterations)
    
    if (contains(version,'2025'))
        queue=createParallelProgressBar_v2025(totalIterations);
    else
        queue=createParallelProgressBar_pre2025(totalIterations);
    end

end