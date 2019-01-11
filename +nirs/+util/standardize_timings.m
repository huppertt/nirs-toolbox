function raw = standardize_timings(raw)
%% This function creates a standardized set of timing info for all files
% Each file will have identical time vector and stimulus events after this is run
% The output stimulus vector is the logical AND of all input stimulus vectors
%
% Files should be aligned beforehand (relative to scan start) typically by
% removing extraneous time using nirs.modules.TrimBaseline
%
% The purpose is to achieve consistency among files, typically for subject-swapping 
% approach to permutation testing for hyperscanning
nsub = length(raw);
conds = sort(raw(1).stimulus.keys);
ncond = length(conds);

%% Determine minimum number of time points across all files (those exceeding will be truncated)
ntpts = inf;
for i = 1:nsub
    ntpts = min(ntpts,length(raw(i).time));
end
time = (0:(ntpts-1))*(1/raw(1).Fs);

%% Find the logical AND of all the stim vectors
X = zeros(ntpts,ncond);
for i = 1:nsub
    tmpX = getTruncatedSortedStims( raw(i) , ntpts , conds );
    X = X + tmpX;
end
X = (X == nsub);

%% Convert common stim vector to stim event
stims = Dictionary();
for i = 1:ncond
    if ~any(X(:,i))
        warning('Empty condition "%s". Check alignment or subjects missing conditions/blocks/trials.',conds{i});
    end
    stims(conds{i}) = nirs.design.vector2event(time,X(:,i),conds{i});
end

%% Set time and stim information equal for everyone
for i=1:nsub
    raw(i).time = time;
    raw(i).stimulus = stims;
    raw(i).data = raw(i).data(1:ntpts,:);
end

end

function X = getTruncatedSortedStims( raw, ntpts, names )

assert(isequal(names,sort(names)),'Names should be sorted');
assert(length(raw)==1,'Only 1 file allowed');

[X,tmpnames] = raw.getStimMatrix;
[tmpnames,sortidx] = sort(tmpnames);
assert(isequal(names,tmpnames),'Stimulus name mismatch');

X = X(1:ntpts,sortidx);

end