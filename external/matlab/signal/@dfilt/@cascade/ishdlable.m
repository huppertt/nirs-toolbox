function [result, errstr, errorObj] = ishdlable(Hb)
%ISHDLABLE True if HDL can be generated for the filter object.
%   ISHDLABLE(Hd) determines if HDL code generation is supported for the
%   filter object Hd and returns true or false.
%
%   The determination is based on the filter structure and the 
%   arithmetic property of the filter.
%
%   The optional second return value is a string that specifies why HDL
%   could not be generated for the filter object Hd.
%
%   See also DFILT, MFILT, GENERATEHDL, GENERATETB.

%   Copyright 2004 The MathWorks, Inc.

  result = true;                        % default setting, can change
  errstr = '';
  errorObj = [];

  first = [];
  for n = 1:length(Hb.Stage)
    % The order of these is important
    if isa(Hb.Stage(n),'dfilt.cascade')
      result = false;
      errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedCascades');
      errstr = getString(errorObj);
      break;
    end
    if ~any(strcmpi(fieldnames(get(Hb.Stage(n))),'Arithmetic'))
      result = false;
      errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedNoArithProps');
      errstr = getString(errorObj);
      break;
    end
    if isempty(first)
      first = Hb.Stage(n).arithmetic;
    end
    if ~strcmpi(Hb.Stage(n).arithmetic, first)
      result = false;
      errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedDiffArithProps');
      errstr = getString(errorObj); 
      break;
    end
    [cando, str, msgObj] = ishdlable(Hb.Stage(n));
    if ~cando
      result = false;
      errstr = str;      
      errorObj = msgObj;
      break;
    end
  end
  %check for farrowsrc in cascade - error when not in last position
  for n = 1:length(Hb.Stage)-1
     if isa(Hb.Stage(n), 'mfilt.farrowsrc')
         result = false;
         errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedFarrowSrc');
         errstr = getString(errorObj);
         break;
     end
  end
  
  if result              % keep testing
      rcf = getratechangefactors(Hb);
      if isa(Hb.Stage(end), 'mfilt.farrowsrc')
          rcf = rcf(1:end-1,:);
      end
      if ~(all(rcf(:,1)==1) || all(rcf(:,2)==1))
          %was any(rcf(:,1)~=1)
          result = false;
          errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedRCFNotMonotonic');
          errstr = getString(errorObj);
      end  
  end
  %keep testing for cascades involving the multrate farrow - only the
  %interp+interp+interp(src) and decim+decim+decim(src) are supported.
  if result && isa(Hb.Stage(end), 'mfilt.farrowsrc')
      isinterp = any(rcf(:,1)~=1);
      orig_rcf = getratechangefactors(Hb);
      if isinterp
          if orig_rcf(end, 1) < orig_rcf(end,2) %
              result = false;
              errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedRCFNotMonotonic');
              errstr = getString(errorObj);
          end
      else % decimating so far
          if orig_rcf(end,1) > orig_rcf(end,2)
              result = false;
              errorObj = message('signal:dfilt:cascade:ishdlable:HdlNotSupportedRCFNotMonotonic');
              errstr = getString(errorObj);
          end
      end
  end
  
% [EOF]

