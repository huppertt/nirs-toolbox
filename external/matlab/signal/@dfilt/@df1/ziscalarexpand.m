function S = ziscalarexpand(Hd,S)
%ZISCALAREXPAND Expand empty or scalar initial conditions to a vector.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

if ~isnumeric(S)
    error(message('signal:dfilt:df1:ziscalarexpand:MustBeNumeric'));
end
if issparse(S),
    error(message('signal:dfilt:df1:ziscalarexpand:Sparse'));
end

numstates = nstates(Hd);

if numstates,
    if isempty(S),
        % If there are states (i.e., numstates > 0) returns a 0
        S = nullstate1(Hd.filterquantizer);
    end

    if length(S)==1 && ~strcmpi(class(S), 'filtstates.dfiir'),
        % Zi expanded to a vector of length equal to the number of states,
        % for both double and single states
        Sinit = S;
        S = filtstates.dfiir;

        % Expand the scalar in Sinit
        S.Numerator   = Sinit(ones(Hd.ncoeffs(1)-1,1));
        S.Denominator = Sinit(ones(Hd.ncoeffs(2)-1,1));
    elseif strcmpi(class(S), 'filtstates.dfiir'),


        if isempty(S.Numerator),
            % Can't use nullstate1 because it returns a FILTSTATES.DFIIR object
            S.Numerator = nullnumstate(Hd.filterquantizer);
        end

        if length(S.Numerator)==1,
            S.Numerator   = S.Numerator(ones(Hd.ncoeffs(1)-1,1));
        end

        if isempty(S.Denominator),
            S.Denominator = nulldenstate(Hd.filterquantizer);
        end
        if length(S.Denominator)==1,
            S.Denominator   = S.Denominator(ones(Hd.ncoeffs(2)-1,1));
        end
        
        lenofs = size(S.Numerator,1) + size(S.Denominator,1);
        if lenofs ~= numstates,
            if size(S.Numerator,1) ~= Hd.ncoeffs(1)-1
                error(message('signal:dfilt:df1:ziscalarexpand:InvalidNumeratorDimensions', Hd.ncoeffs( 1 ) - 1));
            elseif size(S.Denominator,1) ~= Hd.ncoeffs(2)-1
                error(message('signal:dfilt:df1:ziscalarexpand:InvalidDenominatorDimensions', Hd.ncoeffs( 2 ) - 1));
            end
        end

    end

    if ~strcmpi(class(S), 'filtstates.dfiir'),

        % Transpose if row vector only
        if find(size(S)==1),
            S = S(:);
        end

        if size(S,1) ~= numstates,
            error(message('signal:dfilt:df1:ziscalarexpand:InvalidStateDimensions', numstates));
        end

        % Return a FILTSTATES.DFIIR object
        Sinit = S;
        S = filtstates.dfiir;
        S.Numerator = Sinit(1:Hd.ncoeffs(1)-1,:);
        S.Denominator = Sinit(Hd.ncoeffs(1):end,:);
    end
elseif isempty(S) && ~strcmpi(class(S), 'filtstates.dfiir'),
    % Always create a states object (Default filter case).
    Snum = [feval(class(S),zeros(size(S)))];
    Sden = [feval(class(S),zeros(size(S)))];
    S = filtstates.dfiir(Snum,Sden);
else
    if ~isempty(S),
        % This handles the case when there are no states (numstates == 0) and
        % the user attempts to set the states
        error(message('signal:dfilt:df1:ziscalarexpand:MustBeEmpty'));
    end
end
