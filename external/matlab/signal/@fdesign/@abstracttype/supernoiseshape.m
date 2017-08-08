function nsres = supernoiseshape(this,b,linearphase,wl,cb,f,a,args)
%SUPERNOISESHAPE <short description>

%   Copyright 2008-2011 The MathWorks, Inc.


% Process the arguments to set up the respective results structs, get
% the number of iterations to loop, and find which metric to use.
nsres = processargs(b,linearphase,cb,f,a,args);

% Quantize coeffs to WL bits
warnstate = warning('off','fixed:fi:overflow');
bfi = fi(b,1,wl);
fl = bfi.FractionLength;
bq = double(bfi);
bp = bq;
k = nsres.trials;
aux = fi(zeros(k,1),1,wl,fl);

% Initialize the temporary results vector to zeros and the initial
% pareto frontier matrix with the quantized filter results.
bands = noiseshapeparetobands(this);
nbands = length(bands);
m = zeros(k,nbands+1);
p = [zeros(1,nbands) length(b)];
for i=1:nbands,
    % Lowpass, highpass uses {'cb','nb'} while halfband uses {'cb'} for
    % example
    p(i) = measurefreqz(bq,nsres,bands{i}); 
end

% This runs the noise shaping algorithm using a noise shaping
% filter of length 3, then of length 4.  These values have
% experimentally shown to be the best lengths, but future research
% could be done to see if higher values give better results.
refcoeffs   = nsres.refcoeffs;
lsbq = double(lsb(bfi));
nb   = length(b);
n    = nsres.n;
w    = nsres.w;
metric = nsres.metric;
pnorm = nsres.pthnorm;

for nnsf = 3:4
    % The first thing we do is compute a noise shaping filter of
    % the specified length and critical band.  This is generally
    % just a linear phase remez filter, with the pass band set to
    % the critical band.
    nsf = getnoiseshapefilter(this,nnsf,cb);
    
    % This noise shaping algorithm requires that the noise shaping
    % filter be scaled to have the first coefficient be 1, and then
    % have the first coefficient removed from the filter structure.
    nsf = -1*nsf(2:end)'./nsf(1);
    lnsf = length(nsf);
    
    % Initialize the filter state to be random uniform noise
    % within plus or minus 1 LSB.
    r = (rand(k,lnsf)-0.5)*lsbq;
    
    % Initialize the new filter coefficients to zero for now.
    z = zeros(k,n);
    
    % This is the main noise shaper here; essentially, it is
    % just a simple FIR digital filter with a quantization
    % after each time step.  It takes in the original filter
    % coefficient as the input and gives noise shaped filter
    % coefficients as the output.  It uses the 'nsf' (noise
    % shaping filter) to move the quantization noise out of the
    % critical band.
    for j = 1:1:n
        x1 = refcoeffs(j)+r*nsf;
        aux(:) = x1; % This is the fastest way to quantize
        z(:,j) = double(aux);
        r(:,2:end) = r(:,1:(end-1));
        r(:,1) = x1-z(:,j);
    end
    
    % Ok, now that the noise shaping is done we need to compute
    % the metric in both the critical band ('m(:,1)') and the
    % non-critical band(s) ('m(:,2)').
    %
    % These MxN matrices simple compute an M point FFT on an N point
    % input. For example, w1*z computes the M point FFT on N point
    % vector z in the critical band (w2 is for the non-critical band).
    for i = 1:nbands,
         res = nsres.(['res_' bands{i}]);
         m(:,i) = metric(w{i}*z',res.des,pnorm)';
    end
    m(:,end) = length(b);
    
    % If the filter was linear phase, mirror the coefficients.
    if linearphase
        z = [fliplr(z) z(:,(1+mod(nb,2)):end)];
    end
    
    % Only keep the results that lie on the pareto frontier to
    % save memory and time.
    [p,ix] = paretofrontier([p;m]);
    
    % Keep the filter coefficients for points on the pareto
    % front.
    z2 = [bp;z];
    bp = z2(ix,:);
end

nsres = nsresaddresults(nsres,bp,bands);

% Restore warning
warning(warnstate);

end

%--------------------------------------------------------------------------
function nsres = processargs(b,linearphase,cb,f,a,args)

% Define frequency mask
[res1,res2] = definefreqmask(cb,f,a);
grid = [res1.fgrid; res2.fgrid];
des = [res1.des; res2.des];
[grid,ix] = sort(grid,'ascend');
nsres.res.fgrid = grid;
nsres.res.des = des(ix);
nsres.res_cb = res1;
nsres.res_nb = res2;

% Set the number of trials per batch.
nsres.trials = 1000;

% Define metric
switch args.noiseShapeNorm,
    case inf
        nsres.metric = @maximumripple;
    case 2
        nsres.metric = @leastsquares;
    otherwise
        nsres.metric = @leastpthnorm;
end
nsres.pthnorm = args.noiseShapeNorm;

% Compute DFT constants
j = sqrt(-1);
n = size(b,2);
w1 = sum(exp(-1*j*pi*bsxfun(@times,nsres.res_cb.fgrid,0:(n-1))),3);
w2 = sum(exp(-1*j*pi*bsxfun(@times,nsres.res_nb.fgrid,0:(n-1))),3);

if linearphase
    % If the input filter has a linear phase numerator, then the
    % algorithm will preserve linear phase.  The advantage here is
    % that optimizations can be used to make noise shaping run
    % faster.
    odd = mod(n,2);
    n = ceil(n/2);
    refcoeffs = b((n+(~odd)):end);
    if odd
        w1 = w1(:,n:end)+[zeros(size(w1,1),1) fliplr(w1(:,1:(n-1)))];
        w2 = w2(:,n:end)+[zeros(size(w2,1),1) fliplr(w2(:,1:(n-1)))];
    else
        w1 = w1(:,(n+1):end)+fliplr(w1(:,1:n));
        w2 = w2(:,(n+1):end)+fliplr(w2(:,1:n));
    end
else
    % Non-linear phase, so no optimizations.
    refcoeffs = b;
end

nsres.w = {w1;w2};
nsres.refcoeffs = refcoeffs;
nsres.n = n;

end
function [res1,res2] = definefreqmask(cb,f,a)

r = 250;
res.fgrid = zeros(length(f)/2*r,1);
res.des = res.fgrid;

% Mask (assumes piecewise constant amplitude vector a)
for i = 1:length(f)/2
    res.fgrid((((i-1)*r)+1):(i*r)) = linspace(f(2*i-1),f(2*i),r);
    res.des  ((((i-1)*r)+1):(i*r)) = linspace(a(2*i-1),a(2*i),r);
end

% We don't want to measure the frequency response at pi
idx = res.fgrid~=1;
res.fgrid = res.fgrid(idx);
res.des = res.des(idx);
cfgridnum = sum((res.fgrid >= cb(1)) & (res.fgrid <= cb(2)));
res1.fgrid = zeros(cfgridnum,1);
res2.fgrid = zeros(length(res.fgrid)-cfgridnum,1);

fnames = {'fgrid','des'};
for f = 1:1:length(fnames)
    if length(res.(fnames{f})) ~= length(res.fgrid)
        continue;
    end
    
    res1.(fnames{f}) = res1.fgrid;
    res2.(fnames{f}) = res2.fgrid;
    j = 1;
    k = 1;
    for i = 1:1:length(res.fgrid)
        if (res.fgrid(i) >= cb(1)) && (res.fgrid(i) <= cb(2))
            res1.(fnames{f})(j) = res.(fnames{f})(i);
            j = j+1;
        else
            res2.(fnames{f})(k) = res.(fnames{f})(i);
            k = k+1;
        end
    end
end
end
function m = measurefreqz(b,nsres,band)
res = nsres.(['res_' band]);
m = nsres.metric(freqz(b,1,pi*res.fgrid),res.des,nsres.pthnorm);
end
function m = maximumripple(x,des,pnorm) %#ok<INUSD>
% This function determines that maximum difference between array 'x'
% and array 'des'.  It automatically expands singleton dimensions and
% returns an array of values.
m = 20*log10(max(abs(bsxfun(@minus,abs(x),des))))';
end
function m = leastsquares(x,des,pnorm) %#ok<INUSD>
% This function determines the mean square error between array 'x' and
% array 'des'.  It automatically exands singleton dimensions and
% returns an array of values.

m = sqrt(sum(bsxfun(@minus,abs(x),des).^2)'/length(des));
end
function m = leastpthnorm(x,des,pnorm)
% This function determines the P-th norm between array 'x' and array
% 'des'.

L = size(x,2);
m = zeros(L,1);
y = bsxfun(@minus,abs(x),des);
for i = 1:1:L
    m(i) = norm(y(:,i),pnorm);
end
end
function [p,ix] = paretofrontier(d)
% This function computes the pareto frontier of an MxN dataset 'd'.
% It returns the sorted points on the frontier as 'p' and the indexes
% of those points in 'd' as 'ix'.  A point is considered to be on the
% pareto frontier if no other point exists with all N dimensions lower
% than it.

% Compute the pareto frontier efficiently
[d1,m] = unique(d,'rows');
p = true(size(d1,1),1);
for i=1:size(d1,1)
    if p(i)
        p = p & ~(all(bsxfun(@ge,d1,d1(i,:)),2) & ((1:size(d1,1))'~=i));
    end
end
ix = m(p);
p = d1(p,:);
[p1,pix] = sort(p(:,1),'ascend');
p = p(pix,:);
ix = ix(pix);
end
function nsres = nsresaddresults(nsres,bp,bands)
% bp contains filters on the pareto front

% Keep filter with the best attenuation in the critical band (i.e.
% stopband)
nsres.filters.bns = bp(1,:);
nbands = length(bands);
for i = 1:nbands,
    nsres.filters.metrics.(['bns_' bands{i}]) = measurefreqz(bp(1,:),nsres,bands{i}); 
end
nsres.filters = orderfields(nsres.filters);
nsres = orderfields(nsres);
end


% [EOF]
