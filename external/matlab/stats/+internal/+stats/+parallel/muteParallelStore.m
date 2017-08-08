function muteParallelStore( name, value )
%MUTEPARALLELSTORE is a silent STATPARALLELSTORE for use with PCTRUNONALL.
%
%   MUTEPARALLELSTORE is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

    val = internal.stats.parallel.statParallelStore(name, value);
end

