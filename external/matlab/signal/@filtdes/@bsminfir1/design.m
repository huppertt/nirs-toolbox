function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Get frequency specs, they have been prenormalized
Fpass1 = get(d,'Fpass1');
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');
Fpass2 = get(d,'Fpass2');

F = [Fpass1 Fstop1 Fstop2 Fpass2];
A = [1 0 1];
     

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
delta1 = get(d,'Dpass1');
delta2 = get(d,'Dstop');
delta3 = get(d,'Dpass2');
set(d,'magUnits',magUnits);

DEV = [delta1 delta2 delta3];

[N,Wn,BETA,TYPE] = kaiserord(F,A,DEV);

scaleflag = determinescaleflag(d);

b = fir1(N,Wn,TYPE,kaiser(N+1,BETA),scaleflag);

% Construct object
Hd = dfilt.dffir(b);



