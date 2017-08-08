%DESIGNMETHODS  Return the available design methods for a filter designer
%   M = DESIGNMETHODS(D) returns a cell array, M, with the available design
%   methods for the filter designer, D, and its current settings.
%
%   M = DESIGNMETHODS(D, 'default') returns a cell array, M, with the
%   default design method name for the filter designer, D, and its current
%   settings.
%
%   M = DESIGNMETHODS(D, TYPE) returns a cell array, M, with all the FIR
%   design methods available for filter designer, D, when the input TYPE is
%   set to 'fir'. Returns a cell array, M, with all the IIR design methods
%   available for filter designer, D, when the input TYPE is set to 'iir'.
%   By default, when the input TYPE is omitted, all the FIR and IIR design
%   methods are listed.
%
%   M = DESIGNMETHODS(D, 'full') returns a cell array, M, with a list of
%   fully spelled method names available for filter designer, D. For
%   example, instead of 'butter', the list will contain the fully spelled
%   method name: 'Butterworth'.
%
%   M = DESIGNMETHODS(D, ..., 'SystemObject', FLAG) limits the list of
%   design methods to those that have structures that are supported by
%   System objects when FLAG is set to true. M will be empty if no
%   available filter structure is supported by System objects. When the
%   FLAG is set to false, designmethods returns a list of the methods that
%   have filter structures supported by DFILT/MFILT objects. Omitting the
%   'SystemObject', FLAG input pair is equivalent to setting the FLAG to
%   false. (DSP System Toolbox required).
%
%   % EXAMPLE #1 - Construct a lowpass filter designer and check its design methods.
%   d = fdesign.lowpass('N,Fc',10,12000,48000)
%   m = designmethods(d)
%
%   % EXAMPLE #2 - Change the specifications and check the updated methods.
%   d.Specification = 'Fp,Fst,Ap,Ast';
%   m2 = designmethods(d)
%   m3 = designmethods(d, 'iir')
%   m4 = designmethods(d, 'iir', 'full')
%
%   % EXAMPLE #3 - Get help on a particular method.
%   help(d, m2{1})
%
%   See also FDESIGN, FDESIGN/DESIGN, FDESIGN/DESIGNOPTS, FDESIGN/HELP.

%   Copyright 1999-2011 The MathWorks, Inc.


% [EOF]
