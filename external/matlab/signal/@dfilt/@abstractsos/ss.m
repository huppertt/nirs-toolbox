function [A,B,C,D] = ss(Hd)
%SS  Discrete-time filter to state-space conversion.
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 
%
%   See also DFILT.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

warnsv(Hd);

% Build a cascade and invoke the ss method on it
sosM = Hd.sosMatrix;
sv = Hd.ScaleValues;
construct = class(Hd);
construct = construct(1:end-3); % Remove sos from dfilt.df1sos for example
sect = feval(str2func(construct));
sect.Numerator = sosM(1,1:3);
sect.Denominator = sosM(1,4:6);
Hcas = cascade(dfilt.scalar(sv(1)),sect);

% Add sections
N = nsections(Hd);
for n=2:N,
    addstage(Hcas, dfilt.scalar(sv(n)));
    sect = feval(str2func(construct));
    sect.Numerator = sosM(n,1:3);
    sect.Denominator = sosM(n,4:6);
    addstage(Hcas, sect);
end
% addsection(Hcas, dfilt.scalar(sv(N+1)));

[A,B,C,D] = ss(Hcas);

% [EOF]
