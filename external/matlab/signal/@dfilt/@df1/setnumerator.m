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
        warning(message('signal:dfilt:df1:setnumerator:corruptMATFile'));
        
            % The R13 version of the DF1 did not use a circular buffer AND
            % it used to preappend a zero to beginning of the states vector
            % (S(1)) for both the numerator and denoimator.  The R14 version 
            % uses a circular buffer and it also requires that a zero be 
            % prepended to the top of the states for the numerator side
            % only.
            if ~strcmpi(class(this.HiddenStates),'filtstates.dfiir'),
                nb = this.ncoeffs(1);
                na = this.ncoeffs(2);
                Sinit = this.HiddenStates;
                S = filtstates.dfiir;
                
                % Convert linear to circular states (similar to
                % @df1/thissetstates.m)
                tapIndex = this.tapindex;
                if length(tapIndex == 1), tapIndex = [tapIndex 0]; end
                
                % Linear versions of the Num & Den states
                Snum_lin = flipud(Sinit(2:nb,:));
                Snum_lin = [0; Snum_lin];
                Sden_lin = flipud(Sinit(nb+2:end,:));
                Snum = Snum_lin;
                Sden = Sden_lin;
                
                %
                % Numerator states
                %
                zNIdx = tapIndex(1); zNIdx = double(zNIdx);
                Snum(zNIdx+1:nb,:) = Snum_lin(1:nb-zNIdx,:);
                Snum(1:zNIdx,:) = Snum_lin(nb-zNIdx+1:nb,:);

                %
                % Denominator States
                %
                zDIdx = tapIndex(2); zDIdx = double(zDIdx);
                Sden(zDIdx+1:end,:) = Sden_lin(1:na-1-zDIdx,:);
                Sden(1:zDIdx,:) = Sden_lin(na-zDIdx:end,:);

                % Update the FILTSTATES object and store it back in
                % this.Hiddenstates
                S.Numerator = Snum;
                S.Denominator = Sden;
                this.HiddenStates = S;
            end
    end
end

num = dtfwnum_setnumerator(this, num);

% [EOF]
