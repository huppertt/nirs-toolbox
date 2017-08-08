function varargout=sigwin(varargin)
%SIGWIN  Window object.
%
%   SIGWIN window objects are not recommended.
%   Use <a href="matlab:help tocwindows">window functions</a> instead.
%

%   HWIN = SIGWIN.WINDOW(input1,...) returns a window object, HWIN,
%   of type WINDOW. You must specify a window type with SIGWIN. Each
%   window takes one or more inputs. When you specify SIGWIN.WINDOW
%   with no inputs, a window with default property values is created (the
%   defaults depend on the particular window). The properties may
%   be changed with SET(H,PARAMNAME,PARAMVAL).
%
%   SIGWIN.WINDOW can be one of the following:
%
%   <a href="matlab:help sigwin.bartlett">bartlett</a>       - Bartlett window.
%   <a href="matlab:help sigwin.barthannwin">barthannwin</a>    - Modified Bartlett-Hanning window. 
%   <a href="matlab:help sigwin.blackman">blackman</a>       - Blackman window.
%   <a href="matlab:help sigwin.blackmanharris">blackmanharris</a> - Minimum 4-term Blackman-Harris window.
%   <a href="matlab:help sigwin.bohmanwin">bohmanwin</a>      - Bohman window.
%   <a href="matlab:help sigwin.chebwin">chebwin</a>        - Chebyshev window.
%   <a href="matlab:help sigwin.fcloselattopwin">flattopwin</a>     - Flat Top window.
%   <a href="matlab:help sigwin.gausswin">gausswin</a>       - Gaussian window.
%   <a href="matlab:help sigwin.hamming">hamming</a>        - Hamming window.
%   <a href="matlab:help sigwin.hann">hann</a>           - Hann window.
%   <a href="matlab:help sigwin.kaiser">kaiser</a>         - Kaiser window.
%   <a href="matlab:help sigwin.nuttallwin">nuttallwin</a>     - Nuttall defined minimum 4-term Blackman-Harris window.
%   <a href="matlab:help sigwin.parzenwin">parzenwin</a>      - Parzen (de la Valle-Poussin) window.
%   <a href="matlab:help sigwin.rectwin">rectwin</a>        - Rectangular window.
%   <a href="matlab:help sigwin.taylorwin">taylorwin</a>      - Taylor window.
%   <a href="matlab:help sigwin.triang">triang</a>         - Triangular window.
%   <a href="matlab:help sigwin.tukeywin">tukeywin</a>       - Tukey window.
%
% Window object methods.
%   sigwin/generate - Generate a window vector.
%
%   % EXAMPLE: Construct a 128-point chebyshev window
%   Hwin = sigwin.chebwin(128,100);
%   w = generate(Hwin) % Returns the values representing the window 
%   wvtool(Hwin)       % Window Visualization Tool
%
%   For more information, type
%       doc sigwin
%   at the MATLAB command line.
%

%   Copyright 1988-2012 The MathWorks, Inc.


