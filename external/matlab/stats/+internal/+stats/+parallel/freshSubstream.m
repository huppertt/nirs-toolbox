function substream = freshSubstream(S)
%FRESHSUBSTREAM positions S at the start of the next unused substream.
%
%   FRESHSUBSTREAM is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

    stateOnEntry = S.State;
    S.Substream = S.Substream; % this puts S at the start of the Substream
    if ~isequal(S.State,stateOnEntry)
        % S had advanced within the current Substream, so move to the next one
        S.Substream = S.Substream + 1;
    end
    substream = S.Substream;
end
