function Hbest = minimizecoeffwlfir(this,Href,varargin) %#ok<INUSL>
%   This should be a private method

%   Author(s): R. Losada
%   Copyright 2009 The MathWorks, Inc.


fm = getfmethod(Href);

% Initialize
[Hbest,mrfflag] = optimizecoeffwl(Href,varargin{:});

if ~mrfflag,    
    Hd = copy(Hbest); % Copy in case something goes wrong in searchmicoeffwl.
    
    args.Hbest = Hd;
    args.Href = Href;
    try
        Hbest = searchmincoeffwl(fm,args,varargin{:});
    catch ME
        idx = findstr(ME.identifier,':');
        if strcmpi(ME.identifier(idx(end)+1:end),'unsupportedDesignMethod'),
            error(message('signal:dfilt:dtffir:minimizecoeffwlfir:useMaximizeStopband'));
        else
            % Do nothing, return Hbest from above
        end
    end
end
