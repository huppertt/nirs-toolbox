%   The following functions are available for discrete-time filters analysis(type
%   help dfilt/<FUNCTION> to get help on a specific method - e.g. help
%   dfilt/freqz):
%
%   -----------------------------------------------------------------------
%             Analysis 
%   -----------------------------------------------------------------------
%   <a href="matlab:help dfilt/info">info</a>         - Filter information.
%   <a href="matlab:help dfilt/freqz">freqz</a>        - Frequency response of a discrete-time filter.
%   <a href="matlab:help dfilt/phasez">phasez</a>       - Phase response of a discrete-time filter.
%   <a href="matlab:help dfilt/zerophase">zerophase</a>    - Zero-phase response of a discrete-time filter.
%   <a href="matlab:help dfilt/grpdelay">grpdelay</a>     - Group delay of a discrete-time filter.
%   <a href="matlab:help dfilt/phasedelay">phasedelay</a>   - Phase delay of a discrete-time filter.
%   <a href="matlab:help dfilt/impz">impz</a>         - Impulse response of a discrete-time filter.
%   <a href="matlab:help dfilt/impzlength">impzlength</a>   - Length of the impulse response for a discrete-time filter.
%   <a href="matlab:help dfilt/stepz">stepz</a>        - Step response of a discrete-time filter.
%   <a href="matlab:help dfilt/zplane">zplane</a>       - Pole/Zero plot.
%   <a href="matlab:help dfilt/cost">cost</a>         - Cost Estimate. (DSP System Toolbox Required)
%   <a href="matlab:help dfilt/measure">measure</a>      - Measure characteristics of the frequency response. (DSP System Toolbox Required)
%
%   <a href="matlab:help dfilt/order">order</a>        - Filter order.
%   <a href="matlab:help dfilt/coeffs">coeffs</a>       - Filter coefficients.
%   <a href="matlab:help dfilt/fftcoeffs">fftcoeffs</a>    - FFT coefficients (available only with DFILT.FFTFIR)
%   <a href="matlab:help dfilt/firtype">firtype</a>      - Determine the type (1-4) of a linear phase FIR filter.
%   <a href="matlab:help dfilt/tf">tf</a>           - Convert to transfer function.
%   <a href="matlab:help dfilt/zpk">zpk</a>          - Discrete-time filter zero-pole-gain conversion.
%   <a href="matlab:help dfilt/ss">ss</a>           - Discrete-time filter to state-space conversion.
%   <a href="matlab:help dfilt/nsections">nsections</a>    - Number of sections in a discrete-time filter.
%
%   <a href="matlab:help dfilt/isallpass">isallpass</a>    - True for allpass discrete-time filter.
%   <a href="matlab:help dfilt/iscascade">iscascade</a>    - True for cascaded discrete-time filter.
%   <a href="matlab:help dfilt/isfir">isfir</a>        - True for FIR discrete-time filter.
%   <a href="matlab:help dfilt/islinphase">islinphase</a>   - True for linear discrete-time filter.
%   <a href="matlab:help dfilt/ismaxphase">ismaxphase</a>   - True if maximum phase.
%   <a href="matlab:help dfilt/isminphase">isminphase</a>   - True if minimum phase.
%   <a href="matlab:help dfilt/isparallel">isparallel</a>   - True for discrete-time filter with parallel sections.
%   <a href="matlab:help dfilt/isreal">isreal</a>       - True for discrete-time filter with real coefficients.
%   <a href="matlab:help dfilt/isscalar">isscalar</a>     - True if discrete-time filter is scalar.
%   <a href="matlab:help dfilt/issos">issos</a>        - True if discrete-time filter is in second-order sections form.
%   <a href="matlab:help dfilt/isstable">isstable</a>     - True if the filter is stable.
%
%   Second-order sections (DSP System Toolbox Required)
%   <a href="matlab:help dfilt/scale">scale</a>        - Scale second-order sections.
%   <a href="matlab:help dfilt/scalecheck">scalecheck</a>   - Check scale of second-order sections.
%   <a href="matlab:help dfilt/scaleopts">scaleopts</a>    - Create options object for sos scaling.
%   <a href="matlab:help dfilt/specifyall">specifyall</a>   - Fully specify fixed-point filter settings.
%   <a href="matlab:help dfilt/cumsec">cumsec</a>       - Vector of cumulative second-order section filters.
%   <a href="matlab:help dfilt/reorder">reorder</a>      - Reorder second-order sections.
%
%   Fixed-Point (DSP System Toolbox and Fixed-Point Designer Required)
%   <a href="matlab:help dfilt/qreport">qreport</a>      - Quantization report.
%   <a href="matlab:help dfilt/autoscale">autoscale</a>    - Dynamic range scaling.
%   <a href="matlab:help dfilt/set2int">set2int</a>      - Scale the coefficients to integer numbers.
%   <a href="matlab:help dfilt/norm">norm</a>         - Filter norm.
%   <a href="matlab:help dfilt/normalize">normalize</a>    - Normalize filter coefficients.
%   <a href="matlab:help dfilt/denormalize">denormalize</a>  - Restore the coefficients from a NORMALIZE.
%   <a href="matlab:help dfilt/double">double</a>       - Cast filter to double-precision floating-point arithmetic.
%   <a href="matlab:help dfilt/reffilter">reffilter</a>    - Reference double-precision floating-point filter.
%   <a href="matlab:help dfilt/noisepsd">noisepsd</a>     - Power spectral density of filter output due to roundoff noise.
%   <a href="matlab:help dfilt/noisepsdopts">noisepsdopts</a> - Create options object for noisepsd computation.
%   <a href="matlab:help dfilt/freqrespest">freqrespest</a>  - Frequency response estimate via filtering.
%   <a href="matlab:help dfilt/freqrespopts">freqrespopts</a> - Create options object for frequency response estimate.
%
%   Multi-stages 
%   <a href="matlab:help dfilt/cascade">cascade</a>      - Cascade filter objects.
%   <a href="matlab:help dfilt/parallel">parallel</a>     - Connect filters in parallel.
%   <a href="matlab:help dfilt/nstages">nstages</a>      - Number of stages in a cascade or parallel filter.
%   <a href="matlab:help dfilt/addstage">addstage</a>     - Add a stage to a cascade or parallel filter.
%   <a href="matlab:help dfilt/removestage">removestage</a>  - Remove a stage in a cascade or parallel filter.
%   <a href="matlab:help dfilt/setstage">setstage</a>     - Set a stage in a cascade or parallel filter.
%

%   Copyright 2006-2010 The MathWorks, Inc.
