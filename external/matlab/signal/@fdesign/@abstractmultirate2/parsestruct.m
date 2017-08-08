function [filtstruct, varargin] = parsestruct(this, filtstruct, method, varargin)
%PARSESTRUCT   Parse the inputs for the FilterStructure parameter.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.

if nargin > 3
    % Convert structure of options to cell array if necessary
    if isstruct(varargin{1}),
        opts = varargin{1};
        if isfield(opts,'FilterStructure'),
            filtstruct = opts.FilterStructure;
            opts = rmfield(opts,'FilterStructure');
            varargin{1} = opts;
        end
    else
        
        indx = find(strcmpi(varargin, 'filterstructure'));
        if ~isempty(indx)

            % If there are multiple filter structures passed in we just use the
            % last one.
            filtstruct = varargin{indx(end)+1};
            varargin([indx indx+1]) = [];
        end
    end
end

if ~any(strcmpi(filtstruct,validstructures(this,method))),
    error(message('signal:fdesign:abstractmultirate2:parsestruct:invalidStructure'));
end

% [EOF]
