function Y = dbm(varargin)
%DBM Convert to decibels relative to 1 mW (milliWatt).
%   DBM(X) converts the elements of X to decibel units
%   relative to 1mW across a 1 Ohm load.  The elements
%   of X are assumed to represent voltage measurements.
%
%   DBM(X,U) indicates the units of the elements in X,
%   and may be 'power', 'voltage' or any portion of
%   either unit string.  If omitted, U='voltage'.
%
%   DBM(X,R) indicates a measurement reference load of
%   R Ohms.  If omitted, R=1 Ohm.  Note that R is only
%   needed for the conversion of voltage measurements,
%   and is ignored if U is 'power'.
%
%   DBM(X,U,R) specifies both a unit string and a
%   reference load.
%
%   Examples:
%   1) Convert 0.1 volt to dBm (1 Ohm ref.)
%      dbm(.1)           % +10 dBm
%
%   2) Convert sqrt(.5)=0.7071 volts to dBm (50 Ohm ref.)
%      dbm(sqrt(.5),50)  % +10 dBm
%
%   3) Convert 1 mW to dBm
%      dbm(1e-3,'power') % +0 dBm
%
%   See also DB, ABS, ANGLE.

%   Author(s): D. Orofino
%   Copyright 1988-2002 The MathWorks, Inc.

Y=db(varargin{:}) + 30;

% [EOF] dbm.m
