function [coefs,varargout] = robust_cwt(SIG,scales,WAV,plotmode,varargin)
% robust version of the Matlab continious wavelet transform function
% CWT Real or Complex Continuous 1-D wavelet coefficients.
%   COEFS = CWT(S,SCALES,'wname') computes the continuous
%   wavelet coefficients of the vector S at real, positive
%   SCALES, using wavelet whose name is 'wname'.
%   The signal S is real, the wavelet can be real or complex.
%
%   COEFS = CWT(S,SCALES,'wname','plot') computes
%   and, in addition, plots the continuous wavelet
%   transform coefficients.
%
%   COEFS = CWT(S,SCALES,'wname',PLOTMODE) computes and,
%   plots the continuous wavelet transform coefficients.
%   Coefficients are colored using PLOTMODE.
%   PLOTMODE = 'lvl' (By scale) or
%   PLOTMODE = 'glb' (All scales) or
%   PLOTMODE = 'abslvl' or 'lvlabs' (Absolute value and By scale) or
%   PLOTMODE = 'absglb' or 'glbabs' (Absolute value and All scales)
%
%   CWT(...,'plot') is equivalent to CWT(...,'absglb')
%
%   You get 3-D plots (surfaces) using the same keywords listed
%   above for the PLOTMODE parameter, preceded by '3D'. For example:
%   COEFS = CWT(...,'3Dplot') or COEFS = CWT(...,'3Dlvl').
%
%   COEFS = CWT(S,SCALES,'wname',PLOTMODE,XLIM) computes, and
%   plots, the continuous wavelet transform coefficients.
%   Coefficients are colored using PLOTMODE and XLIM.
%   XLIM = [x1 x2] with 1 <= x1 < x2 <= length(S).
%
%   COEFS = CWT(...,'scal') or [COEFS,SC] = CWT(...,'scal') or
%   COEFS = CWT(...,'scalCNT') or [COEFS,SC] = CWT(...,'scalCNT')
%   computes the continuous wavelet transform coefficients and
%   the corresponding scalogram (percentage of energy for each
%   coefficient). When PLOTMODE is equal to 'scal', a scaled
%   image of scalogram is displayed and when PLOTMODE is equal to
%   'scalCNT', a contour representation of scalogram is displayed.
%
%   For each given scale a within the vector SCALES, the
%   wavelet coefficients C(a,b) are computed for b = 1 to
%   ls = length(S), and are stored in COEFS(i,:)
%   if a = SCALES(i).
%   Output argument COEFS is a la-by-ls matrix where la is
%   the length of SCALES. COEFS is a real or complex matrix
%   depending on the wavelet type.
%
%   Examples of valid uses are:
%     t = linspace(-1,1,512);
%     s = 1-abs(t);
%     c = cwt(s,1:32,'cgau4');
%     c = cwt(s,[64 32 16:-2:2],'morl');
%     c = cwt(s,[3 18 12.9 7 1.5],'db2');
%     c = cwt(s,1:32,'sym2','lvl');
%     c = cwt(s,1:64,'sym4','abslvl',[100 400]);
%     [c,Sc] = cwt(s,1:64,'sym4','scal');
%     [c,Sc] = cwt(s,1:64,'sym4','scalCNT');
%
%   See also WAVEDEC, WAVEFUN, WAVEINFO, WCODEMAT.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision 23-May-2012.
%   Copyright 1995-2012 The MathWorks, Inc.

% Check scales.
%--------------
err = 0;
if isempty(scales) ,         err = 1;
elseif min(size(scales))>1 , err = 1;
elseif min(scales)<eps,      err = 1;
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_ScaVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_ScaVal'))
end

% Check signal.
%--------------
if isnumeric(SIG)
    ySIG    = SIG;
    lenSIG  = length(ySIG);
    xSIG    = (1:lenSIG);
    stepSIG = 1;
    
elseif isstruct(SIG)
    try
        ySIG = SIG.y;
    catch ME  %#ok<NASGU>
        err = 1;
    end
    if err~=1
        lenSIG = length(ySIG);
        try
            xSIG = SIG.x; stepSIG = xSIG(2)-xSIG(1);
        catch ME  %#ok<NASGU>
            try
                stepSIG = SIG.step;
                xSIG = (0:stepSIG:(lenSIG-1)*stepSIG);
            catch ME  %#ok<NASGU>
                try
                    xlim = SIG.xlim;
                    xSIG = linspace(xlim(1),xlim(2),lenSIG);
                    stepSIG = xSIG(2)-xSIG(1);
                catch ME  %#ok<NASGU>
                    xSIG = (1:lenSIG); stepSIG = 1;
                end
            end
        end
    end
    
elseif iscell(SIG)
    ySIG = SIG{1};
    xATTRB  = SIG{2};
    lenSIG  = length(ySIG);
    len_xATTRB = length(xATTRB);
    if len_xATTRB==lenSIG
        xSIG = xATTRB;
        stepSIG = xSIG(2)-xSIG(1);
        
    elseif len_xATTRB==2
        xlim = xATTRB;
        xSIG = linspace(xlim(1),xlim(2),lenSIG);
        stepSIG = xSIG(2)-xSIG(1);
        
    elseif len_xATTRB==1
        stepSIG = xATTRB;
        xSIG = (0:stepSIG:(lenSIG-1)*stepSIG);
    else
        xSIG = (1:lenSIG); stepSIG = 1;
    end
else
    err = 1;
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_SigVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_SigVal'))
end

% Check wavelet.
%---------------
getINTEG = 1;
getWTYPE = 1;
if ischar(WAV)
    precis = 10; % precis = 15;
    [val_WAV,xWAV] = intwave(WAV,precis);
    stepWAV = xWAV(2)-xWAV(1);
    wtype = wavemngr('type',WAV);
    if wtype==5 , val_WAV = conj(val_WAV); end
    getINTEG = 0;
    getWTYPE = 0;
    
elseif isnumeric(WAV)
    val_WAV = WAV;
    lenWAV  = length(val_WAV);
    xWAV = linspace(0,1,lenWAV);
    stepWAV = 1/(lenWAV-1);
    
elseif isstruct(WAV)
    try
        val_WAV = WAV.y;
    catch ME  %#ok<NASGU>
        err = 1;
    end
    if err~=1
        lenWAV = length(val_WAV);
        try
            xWAV = WAV.x; stepWAV = xWAV(2)-xWAV(1);
        catch ME  %#ok<NASGU>
            try
                stepWAV = WAV.step;
                xWAV = (0:stepWAV:(lenWAV-1)*stepWAV);
            catch ME  %#ok<NASGU>
                try
                    xlim = WAV.xlim;
                    xWAV = linspace(xlim(1),xlim(2),lenWAV);
                    stepWAV = xWAV(2)-xWAV(1);
                catch ME  %#ok<NASGU>
                    xWAV = (1:lenWAV); stepWAV = 1;
                end
            end
        end
    end
    
elseif iscell(WAV)
    if isnumeric(WAV{1})
        val_WAV = WAV{1};
    elseif ischar(WAV{1})
        precis  = 10;
        val_WAV = intwave(WAV{1},precis);
        wtype = wavemngr('type',WAV{1});
        getINTEG = 0;
        getWTYPE = 0;
    end
    xATTRB  = WAV{2};
    lenWAV  = length(val_WAV);
    len_xATTRB = length(xATTRB);
    if len_xATTRB==lenWAV
        xWAV = xATTRB; stepWAV = xWAV(2)-xWAV(1);
        
    elseif len_xATTRB==2
        xlim = xATTRB;
        xWAV = linspace(xlim(1),xlim(2),lenWAV);
        stepWAV = xWAV(2)-xWAV(1);
        
    elseif len_xATTRB==1
        stepWAV = xATTRB;
        xWAV = (0:stepWAV:(lenWAV-1)*stepWAV);
    else
        xWAV = linspace(0,1,lenWAV);
        stepWAV = 1/(lenWAV-1);
    end
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_WavVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_WavVal'))
end
xWAV = xWAV-xWAV(1);
xMaxWAV = xWAV(end);
if getWTYPE ,  wtype = 4; end
if getINTEG ,  val_WAV = stepWAV*cumsum(val_WAV); end

ySIG   = ySIG(:)';
nb_SCALES = length(scales);
coefs     = zeros(nb_SCALES,lenSIG);
ind  = 1;


err= nirs.math.innovations(ySIG'-mean(ySIG),20);
W =sparse(diag(wfun(err)));
 

for k = 1:nb_SCALES
    a = scales(k);
    a_SIG = a/stepSIG;
    j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));
    if length(j)==1 , j = [1 1]; end
    f            = fliplr(val_WAV(j));
    
    X=-convmtx(f,length(ySIG));
   
    
    yf=ySIG*W*X;
    coefs(ind,:) = -sqrt(a)*wkeep1(diff(yf),lenSIG);
   % coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
    ind          = ind+1;
end

% Test for plots.
%----------------
if nargin<4 , return; end

% Display Continuous Analysis.
%-----------------------------
dummyCoefs = coefs;
NBC = 240;
if strncmpi('3D',plotmode,2)
    dim_plot = '3D';
elseif strncmpi('scal',plotmode,4)
    dim_plot = 'SC';
else
    dim_plot = '2D';
end

if isequal(wtype,5)
    if ~isempty(strfind(plotmode,'lvl'))
        plotmode = 'lvl';
    else
        plotmode = 'glb';
    end
end
switch plotmode
    case {'lvl','3Dlvl'}
        lev_mode  = 'row';   abs_mode  = 0;   msg_Ident = 'By_scale';
        
    case {'glb','3Dglb'}
        lev_mode  = 'mat';   abs_mode  = 0;   msg_Ident = '';
        
    case {'abslvl','lvlabs','3Dabslvl','3Dlvlabs'}
        lev_mode  = 'row';   abs_mode  = 1;   msg_Ident = 'Abs_BS';
        
    case {'absglb','glbabs','plot','2D','3Dabsglb','3Dglbabs','3Dplot','3D'}
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        
    case {'scal','scalCNT'}
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        
    otherwise
        plotmode  = 'absglb';
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        dim_plot  = '2D';
end
if ~isempty(msg_Ident) , msg_Ident = [msg_Ident '_']; end
msg_Ident = [msg_Ident 'Values_of'];

if abs_mode , dummyCoefs = abs(dummyCoefs); end
if nargin==5 && ~isequal(plotmode,'scal') && ~isequal(plotmode,'scalCNT')
    xlim = varargin{1};
    if xlim(2)<xlim(1) , xlim = xlim([2 1]); end
    if xlim(1)<1      , xlim(1) = 1;   end
    if xlim(2)>lenSIG , xlim(2) = lenSIG; end
    indices = xlim(1):xlim(2);
    switch plotmode
        case {'glb','absglb'}
            cmin = min(min(dummyCoefs(:,indices)));
            cmax = max(max(dummyCoefs(:,indices)));
            dummyCoefs(dummyCoefs<cmin) = cmin;
            dummyCoefs(dummyCoefs>cmax) = cmax;
            
        case {'lvl','abslvl'}
            cmin = min(dummyCoefs(:,indices),[],2);
            cmax = max(dummyCoefs(:,indices),[],2);
            for k=1:nb_SCALES
                ind = dummyCoefs(k,:)<cmin(k);
                dummyCoefs(k,ind) = cmin(k);
                ind = dummyCoefs(k,:)>cmax(k);
                dummyCoefs(k,ind) = cmax(k);
            end
    end
elseif isequal(plotmode,'scalCNT')
    if ~isempty(varargin) , nbcl =  varargin{1}; end
end

nb    = min(5,nb_SCALES);
level = '';
for k=1:nb , level = [level ' '  num2str(scales(k))]; end %#ok<AGROW>
if nb<nb_SCALES , level = [level ' ...']; end
nb     = ceil(nb_SCALES/20);
ytics  = 1:nb:nb_SCALES;
tmp    = scales(1:nb:nb*length(ytics));
ylabs  = num2str(tmp(:));
plotPARAMS = {NBC,lev_mode,abs_mode,ytics,ylabs,'',xSIG};

switch dim_plot
    case 'SC'
        if ~exist('nbcl','var') , nbcl = 10; end
        switch plotmode
            case 'scal',     typePLOT = 'image';
            case 'scalCNT' , typePLOT = 'contour';
        end
        SC = wscalogram(typePLOT,coefs,scales,ySIG,xSIG,'nbcl',nbcl);
        if nargout>1 , varargout{1} = SC; end
        
    case '2D'
        if wtype<5
            titleSTR = getWavMSG(['Wavelet:divCMDLRF:' msg_Ident],level);
            plotPARAMS{6} = titleSTR;
            axeAct = gca;
            plotCOEFS(axeAct,dummyCoefs,plotPARAMS);
        else
            axeAct = subplot(2,2,1);
            titleSTR = getWavMSG('Wavelet:divCMDLRF:Real_Part',level);
            plotPARAMS{6} = titleSTR;
            plotCOEFS(axeAct,real(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,2);
            titleSTR = getWavMSG('Wavelet:divCMDLRF:Imag_Part',level);
            plotPARAMS{6} = titleSTR;
            plotCOEFS(axeAct,imag(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,3);
            titleSTR = getWavMSG('Wavelet:divCMDLRF:Modulus_of',level);
            plotPARAMS{6} = titleSTR;
            plotCOEFS(axeAct,abs(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,4);
            titleSTR = getWavMSG('Wavelet:divCMDLRF:Angle_of',level);
            plotPARAMS{6} = titleSTR;
            plotCOEFS(axeAct,angle(dummyCoefs),plotPARAMS);
        end
        colormap(pink(NBC));
        
    case '3D'
        if wtype<5
            titleSTR = getWavMSG(['Wavelet:divCMDLRF:' msg_Ident],level);
            plotPARAMS{6} = titleSTR;
            axeAct = gca;
            surfCOEFS(axeAct,dummyCoefs,plotPARAMS);
        else
            axeAct = subplot(2,2,1);
            titleSTR = ['Real part of Ca,b for a = ' level];
            plotPARAMS{6} = titleSTR;
            surfCOEFS(axeAct,real(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,2);
            titleSTR = ['Imaginary part of Ca,b for a = ' level];
            plotPARAMS{6} = titleSTR;
            surfCOEFS(axeAct,imag(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,3);
            titleSTR = ['Modulus of Ca,b for a = ' level];
            plotPARAMS{6} = titleSTR;
            surfCOEFS(axeAct,abs(dummyCoefs),plotPARAMS);
            axeAct = subplot(2,2,4);
            titleSTR = ['Angle of Ca,b for a = ' level];
            plotPARAMS{6} = titleSTR;
            surfCOEFS(axeAct,angle(dummyCoefs),plotPARAMS);
        end
end



%----------------------------------------------------------------------
function plotCOEFS(axeAct,coefs,plotPARAMS)

[NBC,lev_mode,abs_mode,ytics,ylabs,titleSTR] = deal(plotPARAMS{1:6});

coefs = wcodemat(coefs,NBC,lev_mode,abs_mode);
image(coefs);
set(axeAct, ...
    'YTick',ytics, ...
    'YTickLabel',ylabs, ...
    'YDir','normal', ...
    'Box','On' ...
    );
title(titleSTR,'Parent',axeAct);
xlabel(getWavMSG('Wavelet:divCMDLRF:TimeORSpace'),'Parent',axeAct);
ylabel(getWavMSG('Wavelet:divCMDLRF:Scales_a'),'Parent',axeAct);
%----------------------------------------------------------------------
function surfCOEFS(axeAct,coefs,plotPARAMS)

[NBC,~,~,ytics,ylabs,titleSTR] = deal(plotPARAMS{1:6});

surf(coefs);
set(axeAct, ...
    'YTick',ytics, ...
    'YTickLabel',ylabs, ...
    'YDir','normal', ...
    'Box','On' ...
    );
title(titleSTR,'Parent',axeAct);
xlabel(getWavMSG('Wavelet:divCMDLRF:TimeORSpace'),'Parent',axeAct);
ylabel(getWavMSG('Wavelet:divCMDLRF:Scales_a'),'Parent',axeAct);
zlabel('COEFS','Parent',axeAct);

xl = [1 size(coefs,2)];
yl = [1 size(coefs,1)];
zl = [min(min(coefs)) max(max(coefs))];
set(axeAct,'XLim',xl,'YLim',yl,'ZLim',zl,'view',[-30 40]);

colormap(pink(NBC));
shading('interp')
%----------------------------------------------------------------------


function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
