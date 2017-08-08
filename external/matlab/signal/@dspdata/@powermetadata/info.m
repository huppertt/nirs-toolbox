function varargout = info(this)
%INFO   Return information about the meta data

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

f = {'EstimationMethod', 'FrequencyUnits', 'DataUnits'};
v = {this.FrequencyUnits, lclgetdataunits(this)};

if isempty(this.SourceSpectrum),
    v = {'Unknown', v{:}};
else
    
    % Get all the fields of the source and remove "estimationmethod"
    nf = fieldnames(this.SourceSpectrum);
    if iscell(nf)      
      nf(find(strcmpi(nf,'EstimationMethod'))) = [];
    else          
      nf = {}; 
    end
    f = [f nf(:)'];
    v = {this.SourceSpectrum.EstimationMethod, v{:}};
    
    for indx = 1:length(nf)
        v{end+1} = get(this.SourceSpectrum, nf{indx});
        if isnumeric(v{end}),
            v{end} = num2str(v{end});
        end
    end
end

if nargout > 1,
    varargout = {f, v};
else
    
    for indx = 1:length(f), f{indx} = sprintf('%s:', f{indx}); end
    
    i = [strvcat(f) repmat(' ', length(f), 2) strvcat(v)];

    i = cellstr(i);
    i = sprintf('%s\n', i{:});

    if nargout
        varargout = {i};
    else
        fprintf(1, i);
    end
end

% -------------------------------------------------------------------------
function d = lclgetdataunits(this)

d = this.DataUnits;
if isempty(d),
    d = 'none';
end

% [EOF]
