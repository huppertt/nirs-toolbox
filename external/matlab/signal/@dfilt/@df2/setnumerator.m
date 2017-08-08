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
        warning(message('signal:dfilt:df2:setnumerator:corruptMATFile'));
        
            % The R13 version of the DF2 did not use a circular buffer AND
            % it used to append a zero to the end of the states vector
            % (S(end)).  The R14 version uses a circular buffer and it
            % requires that a zero be prepended to the top of the states
            % vector (S(1)).
            linStates = [0; this.HiddenStates(1:max(ncoeffs)-1)];
            
            % Convert linear to circular states (similar to
            % df2/thissetstates.m)
            tapIndex = this.tapindex(1);
            [nr ncol] = size(linStates);
            Scir = linStates;
            Scir(tapIndex+1:end,:) = linStates(1:nr-tapIndex,:);
            Scir(1:tapIndex,:) = linStates(nr-tapIndex+1:end,:);
            this.HiddenStates = Scir;
    end
end

num = dtfwnum_setnumerator(this, num);

% [EOF]
