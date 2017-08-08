function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.

state = hTar.MapCoeffsToPorts;

if strcmpi(state,'on')
    [mapstate coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);
    % Coefficient names  
    hTar.CoeffNames = coeffnames;
    
    % We need to get the coefficients directly from the filter object
    % rather than from hTar because hTar stores the coefficients from
    % mapcoeffstoports, which does not return the leading coefficient of
    % the denominator (Biquad block does not allow feeding the denominator 
    % through port with the leading coefficients).
    Num = Hd.privNum.';
    Den = Hd.privDen.';
    Den(1,:) = 1./Den(1,:);     % filter structure is 1/a(1)
    g = Hd.ScaleValues;         % this value is always in double.
    
    % create the ScaleValues for fixed-point.
    sv = Hd.privScaleValues;
    if isfi(sv)
        % If the scale values are all one, privScaleValues will return an
        % empty object of fi. We need to restore the scale values back as a
        % fi object when mapcoeffstoports is on.
        g = fi(g,numerictype(sv),fimath(sv));
    end
    
    variables{1} = Num;
    variables{2} = Den;
    variables{3} = g;

    % Coefficient values for export
    setprivcoefficients(hTar,variables);
end

domapcoeffstoports = strcmpi(state,'on');

% [EOF]
