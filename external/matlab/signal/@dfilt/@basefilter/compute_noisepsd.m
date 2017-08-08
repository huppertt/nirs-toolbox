function Hpnn = compute_noisepsd(this,L,varargin)
%COMPUTE_NOISEPSD   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if L ~= round(L) || L < 2,
    error(message('signal:dfilt:basefilter:compute_noisepsd:invalidNTrials'));
end


opts = uddpvparse('dspopts.spectrum', {'noisepsdopts', this}, varargin{:});

M = opts.Nfft;

[Vp,Yp] = nlminputnoutput(this,L,M);
  
sumY2 = sum(Yp.*conj(Yp),2);
sumH = sum(Yp./Vp,2); 
H = sumH/L; % We cannot use freqresp here instead becaue H used in the computation of Pnn below must match Yp,Vp.

% Include abs(Vp) in computation
Pnn = ((sumY2/(L-1)) - max(max(abs(Vp)))^2*real(H.*conj(H))*L/(L-1))/M;

% Compensate for input noise quantization
Pnn = nlminputcomp(get_filterquantizer(this),H,M,Pnn);

Pnn = abs(Pnn); % Make sure it is nonnegative despite roundoff

% Scale Pnn to be a true PSD
args = {'SpectrumType','Twosided','CenterDC',false};
if opts.NormalizedFrequency,
    Pnn = Pnn/(2*pi);    
else
    Pnn = Pnn/opts.Fs;
    args = {args{:},'Fs',opts.Fs};
end

% Construct a psd object. Always make it twosided and nondc-centered.
Hpnn = dspdata.psd(Pnn,args{:});

if ishalfnyqinterval(opts),
    onesided(Hpnn);
elseif opts.CenterDC,
    centerdc(Hpnn);
end


% [EOF]
