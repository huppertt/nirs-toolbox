function Hd = constraincoeffwlfir(this,Href,WL,varargin) %#ok<INUSL>
%CONSTRAINCOEFFWLFIR Constrain coefficient wordlength.
%   This should be a private method

%   Author(s): R. Losada
%   Copyright 2009 The MathWorks, Inc.


fm = getfmethod(Href);

% Try design with noise shaping to see if we meet specs
Hd = optimizecoeffwl(Href,varargin{:});
done = Hd.CoeffWordLength<=WL;


if ~done
    Hbest = Hd;
    
    args.Hbest = Hbest;
    args.Href  = Href;
    args.wl    = WL;
    try
        Hd = searchmincoeffwl(fm,args,varargin{:});
    catch ME
        idx = findstr(ME.identifier,':');
        if strcmpi(ME.identifier(idx(end)+1:end),'unsupportedDesignMethod'),
            error(message('signal:dfilt:dtffir:constraincoeffwlfir:constraincoellwlNotSupported'));
        else
            throwAsCaller(ME);
        end
    end
end

