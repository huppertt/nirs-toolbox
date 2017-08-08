function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.

state = hTar.MapCoeffsToPorts;

if strcmpi(state,'on')
    [mapstate coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);
    
    % create conjugated coefficients when mapcoeffstoports is on
     Lattice = Hd.privlattice.';
    
     % Exported coefficient names and variables
     if isempty(Lattice)
         % This is the header_order0 case for the default object. There will
         % be no coefficient exported to the workspace. Thus, set coeffnames
         % and variables to be empty.
         coeffnames = {};
         variables = {};
     else
         % ceate coefficients when mapcoeffstoports is on
         variables{1} = Lattice;
         % create coefficient name for conjugate coefficients
         coeffnames{2} = [coeffnames{1} 'conj'];
         % create conjugated coefficients when mapcoeffstoports is on
         variables{2} = conj(Lattice);
         coefs = coefficients(Hd);
         if length(coefs{1})==1 && ~usepairinorder0(Hd)
             if useconjugategaininorder0(Hd)
                 variables(1) = [];
                 coeffnames(1) = [];
             else
                 variables(2) = [];
                 coeffnames(2) = [];
             end
         end
     end

    hTar.CoeffNames = coeffnames;
    setprivcoefficients(hTar,variables);
end

domapcoeffstoports = strcmpi(state,'on');

% [EOF]
