function distributeToPool(name,values)
%DISTRIBUTETOPOOL makes an ordered named assignment of values to the parallel pool
%
%   DISTRIBUTETOPOOL(NAME,VALUES) distributes each of the individual values
%   in the cell array VALUES on the client to exactly one of the workers in
%   an open parallel pool.  Subsequently, the variable is accessible by
%   NAME on the worker, using the Statistics and Machine Learning Toolbox
%   function WORKERGETVALUE. The variable value can also be modified on the
%   worker using WORKERUPDATEVALUE.
%
%   The companion function VALUES = RETRIEVEFROMPOOL(NAME) can be used on
%   the client to collect the set of values associated with NAME on each of the 
%   workers.  Together, DISTRIBUTETOPOOL and RETRIEVEFROMPOOL assure matched
%   ordering: for any I, the quantity VALUES{I} returned by RETRIEVEFROMPOOL
%   reflects the quantity VALUES{I} previously supplied to DISTRIBUTETOPOOL,
%   along with any changes to the quantity applied by the worker to which
%   VALUES{I} was distributed.      
%
%   Requirements: DISTRIBUTETOPOOL requires an open parallel pool, and the 
%   length of the cell array VALUES must equal the size of the parallel pool. 
%   Violation of either of these conditions results in an exception.
%
%   See also RETRIEVEFROMPOOL, WORKERGETVALUE, WORKERUPDATEVALUE, STATPARALLELSTORE.  
%
%   DISTRIBUTETOPOOL is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

parfor i=1:internal.stats.parallel.getParallelPoolSize
    % Parfor makes a 1-1 assignment of the elements of streams{}
    % to workers in the parallel pool, but the assignment order is
    % non-deterministic. Keep track of the order that occurred
    % while making the assignment.
    taggedValue = struct('tag',i,'value',values{i});
    internal.stats.parallel.statParallelStore(name,taggedValue);
end


