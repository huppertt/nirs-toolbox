function b = thisgenmcode(h)
%THISGENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs, iargs] = genmcode(h.ResponseTypeSpecs, h);

% Get the constraints from the subclass.
[param, value, desc] = getconstraints(h);
if ~isempty(value),
    params = {params{:}, param};
    values = {values{:}, value};
    descs  = {descs{:},  desc};
    iargs  = sprintf('%s, %s', iargs, param);
end

if isspecify(h),
    
    spfp = convertspecialprops(h);
    if ~isempty(spfp)
        indx = findstr(iargs, ',');
        indx = indx(2);
        iargs = [iargs(1:indx-1) ', fprops'  iargs(indx:end)];
        params = {params{:}, 'fprops'};
        values = {values{:}, genmcodeutils('formatcellstr', spfp)};
        descs  = {descs{:}, 'Special Frequency Properties'};
    end
    
    order = 'N';
    params = {'N', params{:}};
    values = {getmcode(h, 'Order'), values{:}};
    descs  = {'', descs{:}};
else
    
    order = convertorder(h);
    if iscell(order)
        params = {params{:}, 'in'};
        values = {values{:}, getmcode(h, order{2})};
        descs  = {descs{:}, 'Initial order estimate'};
        order  = sprintf('{''%s'', in}', order{1});
    else
        order = sprintf('''%s''', order);
    end
end

params = {params{:}, 'dens'};
values = {values{:}, getmcode(h, 'DensityFactor')};
descs  = {descs{:},  ''};

iargs = sprintf('%s, {dens}', iargs);

% If the Phase is not "unspecified" add it to the parameter list.
phase = get(h, 'Phase');
if ~strcmpi(phase, 'linear')
    params = {params{:}, 'phase'};
    values = {values{:}, sprintf('''%sphase''', lower(phase(1:3)))};
    descs  = {descs{:}, 'Phase Specification'};
    iargs  = sprintf('%s, phase', iargs);
end

% If the FIRType is no "unspecified" add it to the parameter list.
if isdynpropenab(h, 'FIRType'),
    firtype = get(h, 'FIRType');
    if ~strcmpi(firtype, 'unspecified'),
        params = {params{:}, 'firtype'};
        values = {values{:}, sprintf('''%s''', firtype)};
        descs  = {descs{:}, 'FIR Type'};
        iargs  = sprintf('%s, firtype', iargs);
    end
end

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.cr;
b.addcr(designdesc(h));
b.addcr('b  = %s(%s, %s);', designfunction(h), order, iargs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
