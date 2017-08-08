%SCALE  Second-order section scaling.
%   SCALE(Hd) scales the second-order section filter Hd using peak
%   magnitude response scaling (L-infinity), in order to reduce the
%   possibility of overflow when Hd operates in fixed-point arithmetic
%   mode.
%
%   HNEW = SCALE(Hd,Pnorm) If an output is requested a new filter is
%   generated with the scaled second-order sections and the original filter
%   is not changed.
%
%   SCALE(Hd,Pnorm) specifies the norm used to scale the filter. Pnorm can
%   be either a discrete-time-domain norm or a frequency-domain norm. Valid
%   time-domain norms are 'l1','l2', and 'linf'. Valid frequency-domain
%   norms are 'L1','L2', and 'Linf'. Note that L2-norm is equal to l2-norm
%   (Parseval's theorem) but the same is not true for other norms.
%
%   The different norms can be ordered in terms of how stringent they are
%   as follows:
%   'l1' >= 'Linf' >= 'L2' = 'l2' >= 'L1' >= 'linf'.
%   Using the most stringent scaling, 'l1', the filter is the least prone
%   to overflow, but also has the worst signal-to-noise ratio. Linf-scaling
%   is the most commonly used scaling in practice.
%
%   SCALE(Hd,Pnorm,P1,V1,P2,V2,...) specifies optional scaling parameters
%   via parameter-value pairs. Valid pairs are:
%   Parameter               Default     Description/Valid values
%   ---------               -------     ------------------------
%   'sosReorder'            'auto'      Reorder section prior to scaling.
%                                       {'auto','none','up','down','lowpass',
%                                       'highpass','bandpass','bandstop'}
%   'MaxNumerator'          2           Maximum value for numerator
%                                       coefficients
%   'NumeratorConstraint'   'none'      {'none', 'unit', 'normalize','po2'}
%   'OverflowMode'          'wrap'      {'wrap','saturate'}
%   'ScaleValueConstraint'  'unit'      {'unit','none','po2'}
%   'MaxScaleValue'         'Not used'  Maximum value for scale values
%
%   Automatic reordering will only take effect when Hd was obtained as a
%   result from a design using FDESIGN. The sections will be automatically
%   reordered depending on the response type of the design (lowpass,
%   highpass, etc.).
%
%   Note that 'MaxScaleValue' will only be used when 'ScaleValueConstraint'
%   is set to something other than 'unit'. If 'MaxScaleValue' is set to a
%   number, the 'ScaleValueConstraint' will be changed to 'none'.
%
%   If the 'Arithmetic' property of Hd is set to 'double' or 'single', the
%   default values shown will be used for all options that are not
%   specified. If the 'Arithmetic' property is set to 'fixed', the values
%   used for the scaling options will be set according to the settings in
%   Hd. However, if a scaling option is specified that differs from the
%   settings in Hd, this option will be used for scaling purposes but will
%   not change the setting in Hd.
%
%   SCALE(Hd,Pnorm,OPTS) uses an options object to specify the optional
%   scaling parameters in lieu of specifying parameter-value pairs. The
%   OPTS object can be created with OPTS = SCALEOPTS(Hd).
%
%   When the OPTS object input or optional scaling parameter-value pair inputs
%   are omitted, Hd is scaled using the default options specified above.
%
%     % Example: Linf-norm scaling of a lowpass elliptic filter:
%     Hs = fdesign.lowpass; % Create a filter design specifications object.
%     Hd = design(Hs,'ellip');       % Design an elliptic SOS filter
%     scale(Hd,'Linf','ScaleValueConstraint','none','MaxScaleValue',1);
%
%   See also DFILT/SCALEOPTS, DFILT/SCALECHECK, DFILT/REORDER,
%   DFILT/CUMSEC, DFILT/NORM, FDESIGN.

%   Copyright 1988-2009 The MathWorks, Inc.

% [EOF]
