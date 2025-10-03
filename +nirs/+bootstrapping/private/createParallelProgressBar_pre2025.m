function queue = createParallelProgressBar(totalIterations)
    % createParallelProgressBar Initializes a progress bar for parallel
    % computations with dynamic color changing from dark orange to blue.
    %
    % Args:
    %     totalIterations (int): Total number of iterations for the
    %                            progress bar.
    %
    % Returns:
    %     queue (parallel.pool.DataQueue): DataQueue to receive progress
    %                                      updates.
    %
    % Example usage in a parallel loop:
    %     numSamples = 100;
    %     % Create progress bar
    %     queue = createParallelProgressBar(numSamples);
    %     parfor i = 1:numSamples
    %         % Simulate computation
    %         pause(0.1);
    %         % Update progress bar
    %         send(queue, i);
    %     end

    % Initialize DataQueue and Progress Bar
    queue = parallel.pool.DataQueue;
    progressBar = waitbar(0, 'Processing...', 'Name', 'Computation Progress');
    
    % Access the Java-based components of the waitbar
    barChildren = allchild(progressBar);
    
    javaProgressBar = barChildren(1).JavaPeer;  % Access the Java progress bar

    % Enable string painting to show percentage inside the bar
    javaProgressBar.setStringPainted(true);

    % Reset persistent variable count
    persistent count
    count = 0;

    % Define colors between which the bar interpolates
    colorEnd = [12, 123, 220] / 255; % light blue
    colorStart = [171, 94, 0] / 255; % brown-orange

    % Nested function to update progress and color
    function updateProgress(~)
        count = count + 1;
        shareComplete = count / totalIterations;
        
        % Update waitbar position
        waitbar(shareComplete, progressBar);

        % Calculate color transition
        currentColor = (1 - shareComplete) * colorStart + ...
                        shareComplete * colorEnd;
        red = currentColor(1);
        green = currentColor(2);
        blue = currentColor(3);

        % Convert RGB triplet to Java Color
        javaColor = java.awt.Color(red, green, blue);

        % Set the progress bar color
        javaProgressBar.setForeground(javaColor);
        
        % Close progress bar when complete
        if count == totalIterations
            close(progressBar);
            count = [];
        end
    end

    % Add listener to the DataQueue
    afterEach(queue, @updateProgress);
end