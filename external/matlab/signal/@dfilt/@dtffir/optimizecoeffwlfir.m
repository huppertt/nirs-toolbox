function [Hbest,mrfflag] = optimizecoeffwlfir(this,Href,varargin) %#ok<INUSL>
%OPTIMIZECOEFFWL Optimize coefficient wordlength for FIR filters.
%   This should be a private method.

%   Copyright 2009 The MathWorks, Inc.


% Test if Fixed-Point Designer is installed
if ~isfixptinstalled,
     error(message('signal:dfilt:dtffir:optimizecoeffwlfir:fixptTbxRq'));
end

% Parse Inputs
args = wloptiminputparse(Href,varargin{:});
mrfflag= args.MatchRefFilter;

% Minimum wordlength design
try
    minwordlength(Href,args);
catch ME,
    throwAsCaller(ME);
end

Hbest = copy(Href);

if args.noiseShaping
    WL = Hbest.CoeffWordLength;
    try_smaller = true;
    while try_smaller && WL > 2,
        WL = WL-1;
        try_smaller = false;
        try
            % For low-order filters, noise-shaping may change passband gain. The
            % result can be a filter that is not within-spec. We need to add extra
            % checking for such a possibility.
            ng = nominalgain(args.fdesignObj);
            for k = 1:args.NTrials,
                Hns(k) = noiseshape(args.fdesignObj,Href,WL,args); %#ok<AGROW>
                if isspecmet(Hns(k),args.fdesignObj,args) && ...
                        passbandspecmet(args.fdesignObj,Hns(k),ng);
                    
                    try_smaller = true;
                    Hbest = Hns(k);
                    break;
                end
            end
        catch %#ok<CTCH>
            try_smaller = false;
        end
    end
    
    Hbest.privArithmetic = 'fixed'; % Don't set public Arithmetic since that
                                    % will trigger another call to
                                    % minwordlengthfir
    
    % Add one to wordlength since we have decreased by one from the last
    % successful design
    Hbest.CoeffWordLength = WL + 1;
    
    setfdesign(Hbest,getfdesign(Href));  % reset the fdesign and fmethod since
    setfmethod(Hbest, getfmethod(Href)); % they may have been lost during noise
                                         % shaping

end

% [EOF]
