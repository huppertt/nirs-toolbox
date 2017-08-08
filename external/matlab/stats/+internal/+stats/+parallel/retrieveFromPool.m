function values = retrieveFromPool(name)
%RETRIEVEFROMPOOL collects values associated with NAME from parallel pool workers.
%
%   VALUES = RETRIEVEFROMPOOL('NAME') retrieves the value associated with the 
%   string 'NAME' from each of the workers in a parallel pool, and assembles
%   the values on the client in the cell array VALUES.
%   The length of VALUES equals the size of the parallel pool. 
%
%   RETRIEVEFROMPOOL('NAME') must have been preceded by an invocation of
%   DISTRIBUTETOPOOL('NAME',VALUES) with the same parallel pool.
%   Together, DISTRIBUTETOPOOL and RETRIEVEFROMPOOL assure matched ordering:
%   for any I, the quantity VALUES{I} returned by RETRIEVEFROMPOOL
%   reflects the quantity VALUES{I} previously supplied to DISTRIBUTETOPOOL,
%   along with any changes to the quantity applied by the worker to which
%   VALUES{I} was distributed.      
%   
%   Requirements: DISTRIBUTETOPOOL requires an open parallel pool.
%   Violation of this condition results in an exception.
%
%   RETRIEVEFROMPOOL is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.
%
%   See also DISTRIBUTETOPOOL.

%   Copyright 2010-2014 The MathWorks, Inc.

poolsz = internal.stats.parallel.getParallelPoolSize;
permutedValues = cell(poolsz);
values = cell(poolsz);

parfor i=1:poolsz
    permutedValues{i} = internal.stats.parallel.statParallelStore(name);
end

for i=1:poolsz
    % Parfor makes a 1-1 assignment of workers to permutedValues, but the
    % assignment order is non-deterministic. Restore the order that 
    % distributeToPool() recorded.
    tag = permutedValues{i}.tag;
    values{tag} = permutedValues{i}.value;
end



