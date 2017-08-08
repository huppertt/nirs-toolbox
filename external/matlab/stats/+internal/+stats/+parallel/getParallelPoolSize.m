function poolSize = getParallelPoolSize(varargin)
% GETPOOLSIZE: This function returns the size of the parallel pool.
%
% Optional name/value arguments:
%
% Guarded    If TRUE, wrap the call to gcp in a try clause, so as
%            to avoid errors if PCT is unavailable.  If the gcp fails,
%            there is no pool, and we return poolSize = 0;
%            If FALSE (the Default), run gcp naked.  This will error
%            in any context where a PCT license is unavailable.
%
% Nocreate   If TRUE, call gcp('nocreate').  If FALSE (the Default),
%            call gcp without the 'nocreate' option.  In this case, unless 
%            the parallel preferences have been set to inhibit auto-create,
%            this query will open a parallel pool.
%

% Copyright 2013 The MathWorks, Inc.

params =   {'guarded', 'nocreate'};
defaults = { false,     false};
[guarded, nocreate] = internal.stats.parseArgs(params, defaults, varargin{:});

if nocreate
    mygcp = @() gcp('nocreate');
else
    mygcp = @gcp;
end

if guarded
    try
        poolObj = mygcp();
    catch
        poolObj = [];
    end
else
    poolObj = mygcp();
end

if isempty(poolObj)
    poolSize = 0;
else
    poolSize = poolObj.NumWorkers;
end