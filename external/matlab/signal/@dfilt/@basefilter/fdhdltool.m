function varargout = fdhdltool(filterobj, varargin)
%FDHDLTOOL Launch the GUI to generate HDL
%   FDHDLTOOL(Hb) launches the Generate HDL dialog that lets you
%   customize and generate Verilog or VHDL code and test benches for a copy
%   of the quantized filter, Hb. You have to call FDHDLTOOL again if
%   you make changes in the filter, Hb at the command line. 
%
%   The main dialog and two subordinate dialogs provide access to many
%   options for customizing the generated filter and test bench code.
%
%   The Generate HDL main dialog categorizes options in two frames:
% 
%   * HDL Filter - Options you are most likely to customize at least once:
%   language selection, name and location of generated files, reset
%   specifications, and optimizations. 
% 
%   * Test bench types - Options for naming and specifying the type of test
%   bench files to generate, and specific types of stimuli the test benches
%   are to apply.
%
%   EXAMPLE:
%           filtdes = fdesign.lowpass('N,Fc,Ap,Ast',30,0.4,0.05,0.03,'linear');
%           Hb = design(filtdes,'filterstructure','dffir');
%           Hb.arithmetic = 'fixed';
%           fdhdltool(Hb);
%
%   See also GENERATEHDL, GENERATETB, GENERATETBSTIMULUS.

%   Copyright 2006-2012 The MathWorks, Inc.

% check for Filter Design HDL Coder
isfdhdlcinstalled;

[cando, ~, errObj] = ishdlable(filterobj);
if ~cando
    error(errObj);
end

if ~isempty(inputname(1)) % when handle to filterobj is passed
    % make a copy and change the names in the widgets
    filter = filterobj;
    hHdl = fdhdlcoderui.fdhdltooldlg(filter);
    hHdl.setfiltername([inputname(1), '_copy']);
else % when constructor is passed directly
    hHdl = fdhdlcoderui.fdhdltooldlg(filterobj);
    % Use default name
    hHdl.setfiltername(['filter', '_copy']);
end

visState = 'On';
for indx = 1:2:length(varargin)
    if strcmpi(varargin(indx), 'visible')
        visState = varargin{indx+1};
    end
end

if strcmpi(visState, 'On')
    hDlg = DAStudio.Dialog(hHdl);
    if nargout >= 1
        varargout{1} = hDlg;
        if nargout >= 2
            varargout{2} = hHdl;
            if nargout > 2
                warning(message('signal:dfilt:basefilter:fdhdltool:fdhdltool'));
                for i=3:nargout
                    varargout{i} = [];
                end
            end
        end
    else
        % Close the fdhdlcoderui.fdhdltooldlg object if no output required
        l = handle.listener(hHdl, 'CloseDialog', @close_listener);
        set(hHdl, 'HDL_Listener', l);
    end
end

% -------------------------------------------------------------------------
function close_listener(hHDL, en) %#ok
% delete the fdhdlcoderui.fdhdltooldlg object to avoid memory leak.
delete(hHDL);
% clear hHDL;


% [EOF]
