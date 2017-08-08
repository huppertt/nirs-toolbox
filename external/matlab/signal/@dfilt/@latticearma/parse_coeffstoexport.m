function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.

state = hTar.MapCoeffsToPorts;

if strcmpi(state,'on')
    [mapstate coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);

     Lattice = Hd.privlattice.';
    Lattice = Lattice(1:max(find(Lattice~=0)));

    Ladder = Hd.privladder.';
    Ladder = Ladder(1:max(find(Ladder~=0)));

    % create coefficient name for conjugate coefficients
    K = coeffnames{1};
    V = coeffnames{2};
    Kconj = [K 'conj'];
    
    % Zero padding so that the orders of Lattice and Ladder are the same.
    max_order = max(length(Lattice),length(Ladder));
    Lattice = makemaxorder(Lattice,max_order);
    Ladder = makemaxorder(Ladder,max_order);
    
    % Exported coefficient names and variables
    if isempty(Hd.privlattice.')||(max_order==1)
        % When the Lattice is empty or it is the header_order0, we only
        % export the non-conjugated K and V and K.
        coeffnames{1} = K;
        coeffnames{2} = V;
        variables{1} = Lattice;
        variables{2} = Ladder;
    else
        coeffnames{1} = K;
        coeffnames{2} = Kconj;
        coeffnames{3} = V;
        variables{1} = Lattice;
        variables{2} = conj(Lattice);
        variables{3} = Ladder;
    end

    hTar.CoeffNames = coeffnames;
    setprivcoefficients(hTar,variables);
end

domapcoeffstoports = strcmpi(state,'on');

%--------------------------------------------------------------------------
function coeffs = makemaxorder(coeffs,maxorder)
% make the coefficient order equals to MAXORDER by zero padding if needed

currentorder = length(coeffs);
M = abs(maxorder-currentorder);
if currentorder < maxorder
    % If currentorder is less than maxorder, the zero padding is required
    % so that the filter structure is balance.
    coeffs = [coeffs; zeros(M,1)];
elseif currentorder > maxorder
    % If the currentorder (coefficient length) is larger than maxorder,
    % then remove the exceeding coefficients as they are zeros.
    coeffs = coeffs(1:end-M);
end


% [EOF]
