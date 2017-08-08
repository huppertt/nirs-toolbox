function filt2block(b,varargin)
%FILT2BLOCK Generate Simulink filter block
%
%   FILT2BLOCK(B) generates a Discrete FIR Filter block with filter
%   coefficients B.
%
%   FILT2BLOCK(B,'subsystem') generates a Simulink subsystem block that
%   implements FIR filter B using sum, gain, and delay blocks.
%
%   FILT2BLOCK(B,...,'FilterStructure',STRUCT) specify the filter
%   structure, STRUCT, as one of 'directForm' | 'directFormTransposed' |
%   'directFormSymmetric' | 'directFormAntiSymmetric' | 'overlapAdd'. The
%   default is 'directForm'. 'overlapAdd' is only valid when you omit the
%   'subsystem' flag (*).
%
%   FILT2BLOCK(B,A) generates a Discrete Filter block with numerator B and
%   denominator A.
%
%   FILT2BLOCK(B,A,'subsystem') generates a Simulink subsystem block that
%   implements IIR filter numerator B and denominator A using sum, gain,
%   and delay blocks.
%
%   FILT2BLOCK(B,A,...,'FilterStructure',STRUCT) specify the filter
%   structure, STRUCT, as one of 'directForm1' | 'directForm1Transposed' |
%   'directForm2' | 'directForm2Transposed'. The default is 'directForm2'.
%
%   FILT2BLOCK(SOS) generates a Biquad Filter block with second order
%   sections matrix SOS. SOS is a Kx6 matrix, where the number of sections,
%   K, must be greater than or equal to 2. (*) 
%
%   FILT2BLOCK(SOS,'subsystem') generates a Simulink subsystem block that
%   implements biquad filter SOS using sum, gain, and delay blocks.
%
%   FILT2BLOCK(SOS,...,'FilterStructure',STRUCT) specify the filter
%   structure, STRUCT, as one of 'directForm1' | 'directForm1Transposed' |
%   'directForm2' | 'directForm2Transposed'. The default is
%   'directForm2Transposed'.
%
%   FILT2BLOCK(D) generates a filter block with coefficients equal to
%   those available in the digital filter, D. You design a digital filter, 
%   D, by calling the <a href="matlab:help designfilt">designfilt</a> function. 
%   A Discrete FIR Filter block is generated  when D is an FIR filter. A 
%   Biquad Filter block is generated when D is an IIR filter (*).
%
%   FILT2BLOCK(D,'subsystem') generates a Simulink subsystem block that
%   implements the digital filter, D, using sum, gain, and delay blocks.
%
%   FILT2BLOCK(D,...,'FilterStructure',STRUCT) for FIR filters, specify
%   the filter structure, STRUCT, as one of 'directForm' |
%   'directFormTransposed' | 'directFormSymmetric' |
%   'directFormAntiSymmetric' | 'overlapAdd'. The default is
%   'directForm'. 'overlapAdd' is only valid when you omit the
%   'subsystem' flag (*). 
%   For IIR filters, specify the filter structure, STRUCT, as one of
%   'directForm1' | 'directForm1Transposed' | 'directForm2' |
%   'directForm2Transposed'. The default is 'directForm2Transposed'.
%
%   (*) A DSP System Toolbox license is required for this syntax.
%
%   FILT2BLOCK(..., PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) customizes
%   the block with the following options:
%
%   -----------------------  --------------   -----------------------------
%          Parameter              Value                Description
%   -----------------------  --------------   -----------------------------
%   'Destination'            [{'current'}     Specify whether you want to 
%                            'new'            add the block to your current 
%                            <user-defined    Simulink model, create a new 
%                             string>]        model to contain the block, 
%                                             or specify the name of the 
%                                             target subsystem. 
%
%   'BlockName'              <string>         Specify the name for the new 
%                                             block. By default the block
%                                             is named 'filter'.
%
%   'OverwriteBlock'         [true {false}]   Specify if you want to 
%                                             overwrite an existing block 
%                                             with the same name as 
%                                             specified by the 'BlockName' 
%                                             parameter or create a new 
%                                             block.
%
%   'MapCoefficientsToPorts' [true {false}]   Specify if you want to map  
%                                             the filter coefficients to 
%                                             the ports of the block.                                             
%                                            
%   'CoefficientNames'       <cell-array>     Specify the coefficient 
%                                             variable names as a cell
%                                             array of strings.
%                                             'MapCoefficientsToPorts' must 
%                                             be true for this property to
%                                             apply. The default is {'Num'}
%                                             for FIR filters,
%                                             {'Num','Den'} for IIR
%                                             filters, and {'Num','Den','g'} 
%                                             for biquad filters. 
%                                                                                                                 
%   'FrameBasedProcessing'   [{true} false]   Specify if you want the block
%                                             to have frame-based or  
%                                             sample-based input processing. 
%
%   'OptimizeZeros'          [{true} false]   Specify if you want to remove 
%                                             zero-gain blocks. 
%                                             This parameter is only valid 
%                                             when you input the 'subsystem' 
%                                             flag.
% 
%   'OptimizeOnes'           [{true} false]   Specify if you want to replace
%                                             unity-gainblocks with direct 
%                                             connections.
%                                             This parameter is only valid 
%                                             when you input the 'subsystem' 
%                                             flag.
%
%   'OptimizeNegativeOnes'   [{true} false]   Specify if you want to replace 
%                                             negative unity-gain blocks 
%                                             with a sign change at the 
%                                             nearest sum block.
%                                             This parameter is only valid 
%                                             when you input the 'subsystem' 
%                                             flag.
%
%   'OptimizeDelayChains'    [{true} false]   Specify if you want to replace   
%                                             cascaded chains of delay with
%                                             an equivalent single integer
%                                             delay block.
%                                             This parameter is only valid 
%                                             when you input the 'subsystem' 
%                                             flag.
%
%   % Example 1: 
%   %   Generate a Discrete Filter block from IIR coefficients. Set the 
%   %   name of the block to 'Lowpass IIR'.
%
%   [B, A] = butter(6,0.4);
%   filt2block(B,A,'BlockName','Lowpass IIR')
%
%   % Example 2: 
%   %   Generate a Simulink subsystem block that implements an FIR lowpass
%   %   filter using sum, gain, and delay blocks. Specify the input 
%   %   processing to be elements as channels by setting the 
%   %   'FrameBasedProcessing' parameter to false. 
%
%   B = fir1(30,.25);
%   filt2block(B,'subsystem','BlockName','Lowpass FIR',...
%             'FrameBasedProcessing',false)

%   Copyright 2012-2013 The MathWorks, Inc.

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  filtType = 'SOS';
else
  if isempty(varargin) || ischar(varargin{1})
    % Only numerator has been passed, so this is an FIR filter    
    filtType = 'FIR';
  else
    % B,A has been passed so this is a non-SOS IIR filter
    a = varargin{1};    
    if all(size(a)>[1 1])
      error(message('signal:signalanalysisbase:inputnotsupported'));
    end    
    varargin(1) = [];
    filtType = 'IIR';
  end
end

% Look for the 'subsystem' flag, 'Optimize*' inputs, 'on', 'off' flags, 
% boolean inputs, structure inputs, and input processing. 
callBlock = true;
optimizeInputs = false;
filterStructInput = false;
removeInputsIdx = [];

if ~isempty(varargin)
  for idx = 1:length(varargin)
    if ischar(varargin{idx})
      inputVal = lower(varargin{idx});
      switch inputVal
        case{'mapcoeffstoports','optimizenegones','coeffnames','inputprocessing'}
          error(message('signal:signalanalysisbase:invalidParamName',varargin{idx}))   
        case{'optimizezeros','optimizeones','optimizedelaychains'}          
          optimizeInputs = true;
        case 'optimizenegativeones'
          optimizeInputs = true;
          varargin{idx} = 'OptimizeNegOnes';
        case 'mapcoefficientstoports'
          varargin{idx} = 'MapCoeffsToPorts';  
        case 'coefficientnames'
          varargin{idx} = 'CoeffNames';                    
        case 'subsystem'
          callBlock = false;
          removeInputsIdx = [removeInputsIdx idx]; %#ok<*AGROW>
        case 'filterstructure'
          filterStructInput = true;
          removeInputsIdx = [removeInputsIdx idx idx+1];          
          if length(varargin) >= idx + 1
            filterStructureName = varargin{idx+1};
          else
            error(message('signal:signalanalysisbase:unspecifiedFiltStruct'));
          end
        case {'on','off'}
          error(message('signal:signalanalysisbase:invalidOnOffInputs'))   
        case {'framebasedprocessing'}
          if length(varargin) >= idx + 1
            tempValue = varargin{idx+1};
            if ~islogical(tempValue)
              error(message('signal:signalanalysisbase:invalidInputProc'))   
            end
            varargin{idx} = 'InputProcessing';
            if tempValue
              varargin{idx+1} = 'columnsaschannels';
            else
              varargin{idx+1} = 'elementsaschannels';
            end
          else
            error(message('signal:signalanalysisbase:unspecifiedInputProc'));
          end          
      end
    elseif islogical(varargin{idx})
      % Convert logical flags to 'on' or 'off'
      if varargin{idx}
        varargin{idx} = 'on';
      else
        varargin{idx} = 'off';
      end
    end
  end
end

dfiltStructName = [];
if filterStructInput
  dfiltStructName = convertStructName(filterStructureName,filtType);
end

switch filtType
  case 'FIR'
    if isempty(dfiltStructName)
      dfiltStructName = 'dffir';
    end
    if strcmp(dfiltStructName,'fftfir') && ~isfdtbxinstalled
      error(message('signal:signalanalysisbase:invalidFFTFIR'));
    end
    
    Hd = dfilt.(dfiltStructName)(double(b));
    
    if strcmp(dfiltStructName,'fftfir')
      % Make sure we pass a correct block length
      Nfft = Hd.BlockLength+length(b)-1;
      B_Nfft = 2^nextpow2(Nfft);      
      Hd.BlockLength = B_Nfft-length(b)+1;
    end

  case 'IIR'
    if isempty(dfiltStructName)
      dfiltStructName = 'df2';        
    end
    Hd = dfilt.(dfiltStructName)(double(b),double(a));    
  case 'SOS'
    if isempty(dfiltStructName)
      dfiltStructName = 'df2tsos';        
    end
    Hd = dfilt.(dfiltStructName)(double(b));  
end    
    
varargin(removeInputsIdx) = [];

if callBlock
  if strcmp(filtType,'SOS') && ~isfdtbxinstalled
    error(message('signal:signalanalysisbase:noDSTLicenceForSOSBlock'));
  end
  
  if optimizeInputs
    % Optimize inputs are only valid for realizemdl calls
    error(message('signal:signalanalysisbase:invalidOptimizeInputs'));
  end
      
  block(Hd,varargin{:});
  
else
  
  if strcmp(dfiltStructName,'fftfir')
    error(message('signal:signalanalysisbase:invalidFFTFIR'))
  end  
  
  realizemdl(Hd,varargin{:})
end

%--------------------------------------------------------------------------
function dfiltStructName = convertStructName(structName,filtType)

switch lower(structName)
  case 'directform'
    if ~strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'IIR'));
    end
    dfiltStructName = 'dffir';
  case 'directformtransposed'
    if ~strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'IIR'));
    end
    dfiltStructName = 'dffirt';    
  case 'directformsymmetric'
    if ~strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'IIR'));
    end
    dfiltStructName = 'dfsymfir';        
  case 'directformantisymmetric'
    if ~strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'IIR'));
    end
    dfiltStructName = 'dfasymfir';            
  case 'overlapadd'
    if ~strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'IIR'));
    end
    dfiltStructName = 'fftfir';                    
  case 'directform1'
    if strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'FIR'));
    end
    dfiltStructName = 'df1';
    if strcmp(filtType,'SOS')
      dfiltStructName = [dfiltStructName 'sos'];
    end        
  case 'directform1transposed'
    if strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'FIR'));
    end
    dfiltStructName = 'df1t';
    if strcmp(filtType,'SOS')
      dfiltStructName = [dfiltStructName 'sos'];
    end    
  case 'directform2'
    if strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'FIR'));
    end
    dfiltStructName = 'df2';
    if strcmp(filtType,'SOS')
      dfiltStructName = [dfiltStructName 'sos'];
    end    
  case 'directform2transposed'
    if strcmp(filtType,'FIR')
      error(message('signal:signalanalysisbase:invalidFiltStruct',structName,'FIR'));
    end
    dfiltStructName = 'df2t';
    if strcmp(filtType,'SOS')
      dfiltStructName = [dfiltStructName 'sos'];
    end
  otherwise
    error(message('signal:signalanalysisbase:invalidFiltStructAll',structName));
end




