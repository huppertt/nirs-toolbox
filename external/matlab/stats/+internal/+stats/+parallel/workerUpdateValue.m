function workerUpdateValue(name,value)
%WORKERUPDATEVALUE updates the value associated with 'NAME' on a parallel pool worker 
%
%   WORKERUPDATEVALUE('NAME',VALUE) reassigns the value associated with 'NAME'
%   on a parallel pool worker to the new value VALUE.  The preexisting value 
%   will have been initialized on the worker by an invocation of 
%   DISTRIBUTETOPOOL('NAME',VALUES) on the client, and the value may have 
%   subsequently been modified by another invocation of
%   WORKERUPDATEVALUE('NAME',VALUE) on the same worker.
%
%   See also DISTRIBUTETOPOOL, RETRIEVEFROMPOOL, WORKERGETVALUE.  
%
%   WORKERUPDATEVALUE is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

    vstruct = internal.stats.parallel.statParallelStore(name);
    vstruct.value = value;
    internal.stats.parallel.statParallelStore(name,vstruct);
end
