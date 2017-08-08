function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Get frequency specs, they have been prenormalized
Fstop1 = get(d,'Fstop1');
Fpass1 = get(d,'Fpass1');
Fpass2 = get(d,'Fpass2');
Fstop2 = get(d,'Fstop2');

F = [Fstop1 Fpass1 Fpass2 Fstop2];
A = [0 1 0];
     

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
delta1 = get(d,'Dstop1');
delta2 = get(d,'Dpass');
delta3 = get(d,'Dstop2');
set(d,'magUnits',magUnits);

DEV = [delta1 delta2 delta3];

[N,Wn,BETA,TYPE] = kaiserord(F,A,DEV);

scaleflag = determinescaleflag(d);

b = fir1(N,Wn,TYPE,kaiser(N+1,BETA),scaleflag);

% Construct object
Hd = dfilt.dffir(b);



