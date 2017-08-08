function [targs, strs] = convertchoices(this)
%CONVERTCHOICES   Return the structure choices.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Common filter structures for dfilt
strs = {'Direct-Form I',...
    'Direct-Form II',...
    'Direct-Form I Transposed',...
    'Direct-Form II Transposed',...
    'State-Space',...
    'Lattice Autoregressive Moving-Average (ARMA)'};

targs = {'df1','df2','df1t','df2t','statespace','latticearma'};

if issos(this),
    for indx = 1:4
        strs{indx}  = [strs{indx} ', SOS'];
        targs{indx} = [targs{indx} 'sos'];
    end
end

% FIR case
if isfir(this),
    [b,a] = tf(this);
    if a(1)==1,
        strs = {strs{5},...
            'Direct-Form FIR',...
            'Direct-Form FIR Transposed'};
        targs = {targs{5}, 'dffir', 'dffirt'};
    else
        strs = {strs{1:5},...
            'Direct-Form FIR',...
            'Direct-Form FIR Transposed'};
        targs = {targs{1:5}, 'dffir', 'dffirt'};
    end
    if isminphase(this),
        strs = {strs{:},'Lattice Moving-Average Minimum Phase'};
        targs = {targs{:}, 'latticemamin'};
    elseif ismaxphase(this),
        strs = {strs{:},'Lattice Moving-Average Maximum Phase'};
        targs = {targs{:}, 'latticemamax'};
    end
    if islinphase(this),
        ftype = firtype(this);
        switch ftype,
            case {1,2},
                strs = {strs{:},'Direct-Form Symmetric FIR'};
                targs = {targs{:}, 'dfsymfir'};
            case {3,4},
                strs = {strs{:},'Direct-Form Antisymmetric FIR'};
                targs = {targs{:}, 'dfasymfir'};
        end
    end

else % IIR case

    % FD Tbx required for coupled-allpass conversions
    if isfdtbxinstalled,
        strs = {strs{:},...
            'Coupled-Allpass (CA) Lattice',...
            'Coupled-Allpass (CA) Lattice with Power-Complementary (PC) Output'};
        targs = {targs{:}, 'calattice', 'calatticepc'};
    end
    if isallpass(this),
        strs = {strs{:},'Lattice allpass'};
        targs = {targs{:}, 'latticeallpass'};
    end
end

% [EOF]
