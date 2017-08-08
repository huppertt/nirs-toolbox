function Hd = window(this, win, varargin)
%WINDOW Design a window filter.

%   Copyright 1988-2011 The MathWorks, Inc.

d = getdesignobj(this);

if ~isfield(d, 'window')
    error(message('signal:fspecs:hbord:window:invalidDesign'));
end

d = feval(d.window);

% Add the System object property if it applies.
supportedStructs = addsysobjdesignopt(d);

if nargin < 2
    % NO OP
elseif isprop(d, win)
    varargin = {win, varargin{:}}; %#ok<CCAT>
elseif isstruct(win)
    set(d, win);
else

    winstr = getwinstr(win);
    winstr = sprintf('window(h, ''Window'', ''%s'')',winstr);

    warning(message('signal:fspecs:hbord:window:deprecatedFeature', 'Window', winstr, 'help(h, ''window'')'));
    
    % If we have a structure, set it into the object.
    if ~isempty(varargin),
        if isstruct(varargin{1})
            set(d, varargin{1});
            varargin(1) = [];
        end
    end
    set(d, 'Window', win);
end

% If we have any more PV pairs set them into the design object.
if length(varargin) > 1, set(d, varargin{:}); end

% Error out before designing the filter if a System object has been
% requested with a structure that is not supported.
if isprop(d,'SystemObject') && d.SystemObject && ...
    ~any(strcmp(d.FilterStructure,supportedStructs))
    method = 'window';
    error(message('signal:fspecs:basecatalog:SysObjNotSupported',...
      d.FilterStructure,method,'SystemObject',method,'SystemObject'))
end

Hd = design(d, this);
if ishp(this),
    Hd = firlp2hp(Hd);
    % Reset the contained FMETHOD.
    Hd.setfmethod(d);
end


% -------------------------------------------------------------------------
function winstr = getwinstr(win)

if isa(win, 'function_handle')
    winstr = ['@', func2str(win)];
elseif iscell(win)
    winstr = ['{' getwinstr(win{1})];
    for indx = 2:length(win)
        if ischar(win{indx})
            winstr = sprintf('%s, %s', winstr, win{indx});
        else
            winstr = sprintf('%s, %g', winstr, win{indx});
        end
    end
    winstr = [winstr '}'];
elseif isnumeric(win)
    winstr = 'WINDOW_VECTOR';
else
    winstr = ['''' win ''''];
end


% [EOF]
