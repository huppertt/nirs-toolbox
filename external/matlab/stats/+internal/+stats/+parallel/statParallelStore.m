function value = statParallelStore( name, value )
%STATPARALLELSTORE Utility for storing persistent named state.
%   The variable VALUE is stored in a persistent struct with
%   fieldname NAME, if statParallelStore is called with both arguments.
%
%   The variable VALUE is retrieved by calling statParallelStore 
%   with character string NAME as the sole input argument.
%
%   This function can be used to make a variable VALUE resident in
%   the workspaces of parallel workers (aka, labs), and available for 
%   reference in subsequent parfor loops. 
%
%   SPMD Block:
%   If VALUE is a cell array, and statParallelStore() is invoked within 
%   an spmd block, then VALUE{I} will be stored on the lab with labindex I. 
%   If labindex I is greater than the length of the cell array, nothing
%   will be done on labindex I.
%   If VALUE is not a cell array, then VALUE will be stored identically
%   on every lab.
%   
%   PARFOR Loop:
%   If statParallelStore() is called within a parfor loop, then the identical
%   VALUE will be stored on every worker engaged in the loop.
%
%   STATPARALLELSTORE is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.


    persistent STORED_VALUES;
    if nargin > 1
        if iscell(value) && labindex<=length(value)
            STORED_VALUES.(name) = value{labindex};
        else
            STORED_VALUES.(name) = value;
        end
    end
    value = STORED_VALUES.(name);
end

