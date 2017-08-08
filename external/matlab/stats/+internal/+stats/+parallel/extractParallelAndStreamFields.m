function [useParallel, useSubstreams, streams] = extractParallelAndStreamFields( options )
%EXTRACTPARALLELANDSTREAMFIELDS processes parallel computation and stream options. 
%
%   EXTRACTPARALLELANDSTREAMFIELDS is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.


    para = statget(options,'UseParallel');
    sub  = statget(options,'UseSubstreams');
    useParallel     =  strcmpi(para,'always') || (islogical(para)&& para==true);
    useSubstreams   =  strcmpi(sub,'always')  || (islogical(sub) && sub==true);  
    streamArg       =  statget(options,'Streams');
    
    % Repackage the Streams argument
    if isempty(streamArg)
        streams = {};
        if useSubstreams
            % If useSubstreams is true, reproducibility requires a single
            % stream to be used within any loop, even if many workers execute
            % iterations of the loop in parallel.  If no specific stream was
            % supplied in the 'options' argument, we take a snapshot
            % of the current default stream on the client, and use this as 
            % the basis for random number generation, even if the default
            % stream subsequently changes.
            streams{1} = RandStream.getGlobalStream;
        end
    elseif ~iscell(streamArg)   % we handle stream arguments with a cell array
        streams = {streamArg};
    else
        streams = streamArg;
    end
    
end %-extractParallelAndStreamFields
