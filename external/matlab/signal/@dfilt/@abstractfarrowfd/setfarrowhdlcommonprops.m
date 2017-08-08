function setfarrowhdlcommonprops(this, hF)
%SETFARROWHDLCOMMONPROPS Set the farrowhdlcommonprops
%   OUT = SETFARROWHDLCOMMONPROPS(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

arith = this.Arithmetic;
hF.arithmetic = arith;

coeffs = coefficients(this);
hF.Coefficients = coeffs{1};

arithisdouble = strcmpi(arith, 'double');
if arithisdouble
    [ivtype, isltype] = conv2hdlsharedtypes(arith);
    [fdvtype, fdsltype] = conv2hdlsharedtypes(arith);    
    [ovtype, osltype] = conv2hdlsharedtypes(arith);
    [iportvtype, iportsltype] = conv2hdlsharedtypes(arith);
    [oportvtype, oportsltype] = conv2hdlsharedtypes(arith);   
    [fdportvtype, fdportsltype] = conv2hdlsharedtypes(arith);     
    % input and output ports need to be adjusted for the vtypes and sltypes
    if hdlgetparameter('isverilog')
        iportvtype = 'wire [63:0]';               % ports are wires
        fdportvtype = 'wire [63:0]';
        oportvtype = 'wire [63:0]';   
    else
        iportvtype = 'real';
        fdportvtype = 'real';
        oportvtype = 'real';        
    end 
    [pvtype, psltype] = conv2hdlsharedtypes(arith);
    [avtype, asltype] = conv2hdlsharedtypes(arith);    
    rmode = 'floor';
    omode = false;
else
    % all the quantizer are *always* signed except for coefficients that
    % can be unsigned too. % fd is always "unsigned"
    
    [iportvtype, iportsltype] = conv2hdlsharedtypes(arith, this.InputWordlength, this.InputFraclength, true);
    [fdportvtype, fdportsltype] = conv2hdlsharedtypes(arith, this.FDWordlength, this.FDFraclength, false);    
    % input and output ports need to be adjusted for the vtypes and sltypes
    if hdlgetparameter('filter_input_type_std_logic') == 1
        [iportvtype, iportsltype] = hdlgetporttypesfromsizes(this.InputWordlength, this.InputFraclength, true);
        [fdportvtype, fdportsltype] = hdlgetporttypesfromsizes(this.FDWordlength, this.FDFraclength, false);
    end

    [oportvtype, oportsltype] = conv2hdlsharedtypes(arith, this.OutputWordlength, this.OutputFraclength, true);
    if hdlgetparameter('filter_output_type_std_logic') == 1
        [oportvtype, oportsltype] = hdlgetporttypesfromsizes( this.OutputWordlength, this.OutputFraclength, true);
    end
    [ivtype, isltype] = conv2hdlsharedtypes(arith, this.InputWordlength, this.InputFraclength, true);
    [fdvtype, fdsltype] = conv2hdlsharedtypes(arith, this.FDWordlength, this.FDFraclength, false);    
    [ovtype, osltype] = conv2hdlsharedtypes(arith, this.OutputWordlength, this.OutputFraclength, true);
    [pvtype, psltype] = conv2hdlsharedtypes(arith, this.ProductWordlength, this.ProductFraclength, true);
    [avtype, asltype] = conv2hdlsharedtypes(arith, this.AccumWordlength, this.AccumFraclength, true);
    rmode = this.Roundmode;
    omode = strcmpi(this.overflowmode, 'saturate');
end

hF.inputportvtype = iportvtype;
hF.inputportsltype = iportsltype;

hF.FDportvtype = fdportvtype;
hF.FDportsltype = fdportsltype;

hF.Outputportvtype = oportvtype;
hF.Outputportsltype = oportsltype;

hF.inputvtype = ivtype;
hF.inputsltype = isltype;

hF.FDvtype = fdvtype;
hF.FDsltype = fdsltype;

hF.Outputvtype = ovtype;
hF.Outputsltype = osltype;

hF.productvtype = pvtype;
hF.productsltype = psltype;

hF.accumvtype = avtype;
hF.accumsltype = asltype;

hF.Roundmode = rmode;
hF.FilterStructure = this.FilterStructure;
hF.Overflowmode = omode;

hF.Comment = info(this);

% [EOF]
