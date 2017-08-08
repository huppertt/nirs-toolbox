function [hTar,domapcoeffstoports] = parse_coeffstoexport(Hd,hTar)
%PARSE_COEFFSTOEXPORT Store coefficient names and values into hTar for
%export.

%   Copyright 2009 The MathWorks, Inc.

% check if spblks is needed for overall multistage filter (dfilt or mfilt)
checkifspblksisneeded(Hd);

userdefinedflag = false;
state = hTar.MapCoeffsToPorts;

if ~isempty(hTar.CoeffNames)
    userdefinedflag = true;
end

% check if all sections is realizable
nsections = length(Hd.Stage);
for k=1:nsections,
    if ~isrealizable(Hd.Stage(k))
        error(message('signal:dfilt:multistage:parse_coeffstoexport:Notsupported', Hd.Stage( k ).FilterStructure));
    end
end
    

if strcmpi(state,'on')
    [~, coeffnames var] = mapcoeffstoports(Hd,'MapCoeffsToPorts','on',...
                                        'CoeffNames',hTar.CoeffNames);
    hTar.CoeffNames = coeffnames;
    setprivcoefficients(hTar,var);
end

domapcoeffstoports = strcmpi(state,'on');

% Append stage index to coefficient names if the coefficient names are not
% user-defined to prevent repeated names.
if ~userdefinedflag && domapcoeffstoports
    coeffnames = hTar.CoeffNames;
    fd = fields(coeffnames);
    for stage = 1:length(fd)
        stagecoeff = coeffnames.(sprintf('Stage%d',stage));
        coeffname_temp = appendcoeffstageindex(Hd,stagecoeff,num2str(stage));    % add stage number
        coeffnames.(sprintf('Stage%d',stage)) = coeffname_temp;
    end
    hTar.CoeffNames = coeffnames;
end



% [EOF]
