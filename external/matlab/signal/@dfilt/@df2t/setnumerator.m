function num = setnumerator(this, num)
%SETNUMERATOR   Set the numerator.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

oldlength = 0;
ncoeffs   = this.ncoeffs;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end

if oldlength~=length(this.Numerator)
    % We are in a load state. Check to see if we've saved a bad object.
    if length(this.HiddenStates) >= max(ncoeffs)
        warning(message('signal:dfilt:df2t:setnumerator:corruptMATFile'));
        this.HiddenStates = this.HiddenStates(1:max(ncoeffs)-1);
    end
end

num = dtfwnum_setnumerator(this, num);

% [EOF]
