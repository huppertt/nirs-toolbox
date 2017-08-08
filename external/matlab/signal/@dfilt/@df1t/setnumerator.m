function num = setnumerator(this, num)
%SETNUMERATOR   Set the numerator.

%   Author(s): P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

oldlength = 0;
ncoeffs   = this.ncoeffs;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end

if oldlength~=length(this.Numerator)
    % We are in a load state. Check to see if we've saved a bad object.
    if length(this.HiddenStates) >= max(ncoeffs)
        warning(message('signal:dfilt:df1t:setnumerator:corruptMATFile'));
        
            % The R13 version of the DF1T used to append a zero to end of 
            % the states vector (S(end)) for both the numerator and denoimator.  
            % The R14 version does not require any zero-padding of the
            % states.
            if ~strcmpi(class(this.HiddenStates),'filtstates.dfiir'),
                nb = this.ncoeffs(1);
                na = this.ncoeffs(2);
                Sinit = this.HiddenStates;
                S = filtstates.dfiir;
                % Negate the denominator states to ensure that output
                % signal matches the R13 result
                S.Denominator = -(Sinit(1:na-1,:));
                S.Numerator = Sinit(na+1:end-1,:);
                this.HiddenStates = S;
            end
    end
end

num = dtfwnum_setnumerator(this, num);

% [EOF]
