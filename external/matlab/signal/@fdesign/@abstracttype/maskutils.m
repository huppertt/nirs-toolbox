function fcns = maskutils(~, isconstrained, axes, units, fs, freqscale, xlim)
%MASKUTILS  Utilities for drawing the masks.
%   MASKUTILS(H, ISC, HAXES, UNITS, FS)  Returns a structure of function
%   handles to help draw masks.  

%   Copyright 2005-2012 The MathWorks, Inc.

% This should be private

narginchk(2,7);

% Parse the inputs
if nargin < 7
  xlim = [];
  if nargin < 6
    freqscale = 'linear';
    if nargin < 5
      fs = 2;
      if nargin < 4
        units = 'db';
        if nargin < 3
          axes = -1;
        end
      end
    end
  end
end

if ischar(axes)
    if nargin > 3
        fs = units;
    end
    units = axes;
    axes = -1;
end

% Return a structure of function handles so that the caller has access to
% this functionality with the passed settings.
fcns.formatapass  = @formatapass;
fcns.formatastop  = @formatastop;
fcns.gethighlow   = @gethighlow;
fcns.getarbmag    = @getarbmag;
fcns.getarbphase  = @getarbphase;
fcns.findbottom   = @findbottom;
fcns.getfs        = @getfs;
fcns.getunits     = @getunits;
fcns.getfreqscale = @getfreqscale;
fcns.getxlim      = @getxlim;
% -------------------------------------------------------------------------
    function out_fs = getfs
        out_fs = fs;
        if isempty(out_fs)
          out_fs = 2;
        end
    end
% -------------------------------------------------------------------------
    function out_units = getunits
        out_units = units;
    end
% -------------------------------------------------------------------------
    function out_freqscale = getfreqscale
        out_freqscale = freqscale;
    end
% -------------------------------------------------------------------------
    function out_xlim = getxlim
        out_xlim = xlim;
    end  
% -------------------------------------------------------------------------
    function A = gethighlow(specs)
        % Returns the magnitude information for the highpass and lowpass
        % case.  They will be differentiated by their frequency edges.

        % Get the attenuation and ripples into the correct units.
        apass = formatapass(specs.Apass);
        astop = formatastop(specs.Astop);
        
        A = [astop(1:2) apass(1) apass(1) apass(2) apass(2) astop(3:4)];
        if apass(1) == apass(2)
            A(4) = nan;
        end
    end
% -------------------------------------------------------------------------
    function A = getarbmag(A)
        
          A = convertmagunits(A, 'linear', lower(units), 'amplitude');
    end
% -------------------------------------------------------------------------
    function P = getarbphase(P)
        
        if strcmpi(units,'degrees'),
          P = convert2deg(P);
        end
    end
% -------------------------------------------------------------------------
    function astop = formatastop(astop)
        % Returns the Astop information in a 4 element vector.
        % [Astop(top) Astop(top)(at edge) Astop(bottom)(at edge)
        % Astop(bottom)]

        if isdb
            
            % If we are not given an Astop, we want to draw the Fstop line
            % down to the bottom of the plot when we are in dB.
            if isempty(astop) || isnan(astop)
                astop = [NaN repmat(findbottom(nan), 1, 2) NaN];
            else
                astop = [-astop -astop findbottom(-astop) NaN];
            end
        else
            
            % If we are not given an Astop, we want to draw the Fstop line
            % to 0 when not in dB.
            if isempty(astop) || isnan(astop)
                astop = [0 0 0 0];
            else
                astop = convertmagunits(astop, 'db', getmagunits, 'stop');
                if strcmpi(units, 'zerophase')
                    astop = [astop astop -astop -astop];
                else
                    astop = [astop astop 0 0];
                end
            end
        end
    end
% -------------------------------------------------------------------------
    function apass = formatapass(apass)

        % Return both the top and the bottom of the box.  This will take
        % into account the 'isconstrained' setting.
        if isdb
            if isempty(apass) || isnan(apass)
                apass = [0 0];
            elseif isconstrained
                apass = [0 -apass];
            else
                apass = [apass -apass]/2;
            end
        else
            if isempty(apass) || isnan(apass)
                apass = [1 1];
            else
                apass = convertmagunits(apass, 'db', getmagunits, 'pass');
                if isconstrained
                    if strcmpi(units, 'squared')
                        apass = [1 apass];
                    else
                        apass = [1 1-apass*2];
                    end
                else
                    if strcmpi(units, 'squared')
                        apass = (1-apass)/2;
                    end
                    apass = 1+[apass -apass];
                end
            end
        end
    end

% -------------------------------------------------------------------------
    function b = isdb
        
        b = strcmpi(units, 'db');
        
    end

% -------------------------------------------------------------------------
    function magunits = getmagunits
        
        magunits = lower(units);

        if strcmpi(magunits, 'zerophase')
            magunits = 'linear';
        end
        
    end

% -------------------------------------------------------------------------
    function b = findbottom(minvalue)

        % If the axes is not a handle (has not been passed) we use a bottom
        % of -60 dB or 0 for linear/squared/zerophase.
        if ishandle(axes)
        
            ylim = get(axes, 'ylim');
            b    = ylim(1);
        elseif isdb
            b = -60;
        else
            b = 0;
        end

        if ~isnan(minvalue)
            % Give it a buffer of 20% so if they resize with matlab autoscaling
            % they won't see the bottom.
            b = min(b, minvalue)*1.2;
        end
    end
end

% [EOF]
