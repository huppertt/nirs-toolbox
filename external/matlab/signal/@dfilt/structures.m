%   DFILT.<STRUCTURE> can be one of the following (type help dfilt.<STRUCTURE>
%   to get help on a specific structure - e.g. help dfilt.df1):
%
%   Note that one usually does not construct DFILT filters explicitly.
%   Instead, one obtains these filters as a result from a design using
%   FDESIGN. 
%     
%   -----------------------------------------------------------------------
%             FIR 
%   -----------------------------------------------------------------------
%   <a href="matlab:help dfilt.dffir">dffir</a>             - Direct-form FIR
%   <a href="matlab:help dfilt.dffirt">dffirt</a>            - Direct-form FIR transposed
%   <a href="matlab:help dfilt.dfsymfir">dfsymfir</a>          - Direct-form symmetric FIR
%   <a href="matlab:help dfilt.dfasymfir">dfasymfir</a>         - Direct-form antisymmetric FIR
%   <a href="matlab:help dfilt.fftfir">fftfir</a>            - Overlap-add FIR
%   <a href="matlab:help dfilt.latticemamax">latticemamax</a>      - Lattice moving-average (MA) for maximum phase
%   <a href="matlab:help dfilt.latticemamin">latticemamin</a>      - Lattice moving-average (MA) for minimum phase
%   <a href="matlab:help dfilt.farrowlinearfd">farrowlinearfd</a>    - Farrow linear fractional delay (*)
%   <a href="matlab:help dfilt.farrowfd">farrowfd</a>          - Farrow fractional delay (*)
%
%   -----------------------------------------------------------------------
%             IIR 
%   -----------------------------------------------------------------------
%   <a href="matlab:help dfilt.allpass">allpass</a>           - Minimum-multiplier allpass filter (*)
%   <a href="matlab:help dfilt.wdfallpass">wdfallpass</a>        - Wave digital allpass filter (*)
%   <a href="matlab:help dfilt.df1">df1</a>               - Direct-form I
%   <a href="matlab:help dfilt.df1sos">df1sos</a>            - Direct-form I, second-order sections
%   <a href="matlab:help dfilt.df1t">df1t</a>              - Direct-form I transposed
%   <a href="matlab:help dfilt.df1tsos">df1tsos</a>           - Direct-form I transposed, second-order sections
%   <a href="matlab:help dfilt.df2">df2</a>               - Direct-form II
%   <a href="matlab:help dfilt.df2sos">df2sos</a>            - Direct-form II, second-order sections
%   <a href="matlab:help dfilt.df2t">df2t</a>              - Direct-form II transposed
%   <a href="matlab:help dfilt.df2tsos">df2tsos</a>           - Direct-form II transposed, second-order sections
%   <a href="matlab:help dfilt.latticeallpass">latticeallpass</a>    - Lattice allpass
%   <a href="matlab:help dfilt.latticear">latticear</a>         - Lattice autoregressive (AR)
%   <a href="matlab:help dfilt.latticearma">latticearma</a>       - Lattice autoregressive moving-average (ARMA)
%   <a href="matlab:help dfilt.statespace">statespace</a>        - State-space
%
%   -----------------------------------------------------------------------
%             Multi-stages 
%   -----------------------------------------------------------------------
%   <a href="matlab:help dfilt.delay">delay</a>             - Integer delay
%   <a href="matlab:help dfilt.scalar">scalar</a>            - Scalar
%   <a href="matlab:help dfilt.cascade">cascade</a>           - Cascade (filters arranged in series)
%   <a href="matlab:help dfilt.parallel">parallel</a>          - Parallel (filters arranged in parallel)
%   <a href="matlab:help dfilt.cascadeallpass">cascadeallpass</a>    - Cascade of minimum-multiplier allpass filters (*)
%   <a href="matlab:help dfilt.cascadewdfallpass">cascadewdfallpass</a> - Cascade of wave digital allpass filters (*)
%   <a href="matlab:help dfilt.calattice">calattice</a>         - Coupled-allpass (CA) lattice (*)
%   <a href="matlab:help dfilt.calatticepc">calatticepc</a>       - Coupled-allpass (CA) lattice with power complementary (PC) output (*)                           
%
%   (*) DSP System Toolbox required
%
%   See also DFILT

%   Copyright 2005-2010 The MathWorks, Inc.



% [EOF]
