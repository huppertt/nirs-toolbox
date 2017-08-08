function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.

state = hTar.MapCoeffsToPorts;

if strcmpi(state,'on')
    [mapstate coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);

    % Coefficient values for export
    num = var{1};
    den = var{2}; den(1) = 1/den(1); % filter structure is 1/a(1)
    
    variables{1} = num;
    variables{2} = den;
    
    % Coefficient names and variables                            
    hTar.CoeffNames = coeffnames;
    setprivcoefficients(hTar,variables);

end

domapcoeffstoports = strcmpi(state,'on');




% [EOF]
