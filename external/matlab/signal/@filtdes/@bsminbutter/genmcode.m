function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

b = sigcodegen.mcodebuffer;

[p,v] = abstract_genmcode(h, d);
p{end+1} = 'match';
v{end+1} = sprintf('''%s''', get(d, 'MatchExactly'));

b.addcr(b.formatparams(p,v));
b.cr;
b.addcr(designdesc(d));
b.addcr('h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2%s);', getfsinput(d));
b.add('Hd = design(h, ''butter'', ''MatchExactly'', match);');

% [EOF]
