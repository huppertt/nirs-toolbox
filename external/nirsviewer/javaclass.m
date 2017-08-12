% Return java.lang.Class instance for MatLab type.
%
% Input arguments:
% mtype:
%    the MatLab name of the type for which to return the java.lang.Class
%    instance
% ndims:
%    the number of dimensions the MatLab data type
% method:
%    the conversion method, "Ragged" to use Java-style multidimensional
%    (ragged) arrays (i.e. arrays of arrays), "Rectangular" to use a
%    more efficient Java implementation where data is stored contiguously
%    (as in MatLab)

% Copyright 2009 Levente Hunyadi
function jclass = javaclass(mtype, ndims, method)

validateattributes(mtype, {'char'}, {'nonempty','row'});
if nargin < 2
    ndims = 0;
else
    validateattributes(ndims, {'numeric'}, {'nonnegative','integer','scalar'});
end
if nargin < 3
    method = 'Ragged';
else
    method = validatestring(method, {'Ragged','Rectangular'});
end

if ndims > 0
    switch method
        case 'Ragged'
            jclassname = javaarrayclass(mtype, ndims);
        case 'Rectangular'
            switch ndims
                case 2
                    jclassname = hu.bme.aut.matlab.Matrix.class;
                case 1
                    jclassname = javaarrayclass(mtype, 1);
                otherwise
                    warning('java:javaclass:NotImplemented', 'Multidimensional array support is not implemented yet.');
                    jclassname = hu.bme.aut.matlab.NumericArray.class;
            end
    end
else
    % The static property .class applied to a Java type returns a string in
    % MatLab rather than an instance of java.lang.Class. For this reason,
    % use a string and java.lang.Class.forName to instantiate a
    % java.lang.Class object; the syntax java.lang.Boolean.class will not
    % do.
    switch mtype
        case 'logical'  % logical array of true and false values
            jclassname = 'java.lang.Boolean';
        case 'char'  % character array
            jclassname = 'java.lang.Character';
        case {'int8','uint8'}  % 8-bit signed and unsigned integer array
            jclassname = 'java.lang.Byte';
        case {'int16','uint16'}  % 16-bit signed and unsigned integer array
            jclassname = 'java.lang.Short';
        case {'int32','uint32'}  % 32-bit signed and unsigned integer array
            jclassname = 'java.lang.Integer';
        case {'int64','uint64'}  % 64-bit signed and unsigned integer array
            jclassname = 'java.lang.Long';
        case 'single'  % single-precision floating-point number array
            jclassname = 'java.lang.Single';
        case 'double'  % double-precision floating-point number array
            jclassname = 'java.lang.Double';
        case 'complexsingle'
            jclassname = 'hu.bme.aut.matlab.ComplexF';
        case 'complexdouble'
            jclassname = 'hu.bme.aut.matlab.Complex';
        otherwise
            error('java:javaclass:InvalidArgumentValue', ...
                'MatLab type "%s" is not supported in Java.', mtype);
    end
end
% When querying a java.lang.Class object by name with the method
% jclass = java.lang.Class.forName(jclassname);
% MatLab generates an error. On one hand, MatLab requires class loader to
% be specified explicitly for the Class.forName method to work. On the
% other hand, the class has to be present on the static class path, classes
% on the dynamic class path will be ignored. Edit classpath.txt and add the
% class to the static class path.
jclass = java.lang.Class.forName(jclassname, true, java.lang.Thread.currentThread().getContextClassLoader());

function jclassname = javaarrayclass(mtype, ndims)

switch mtype
    case 'logical'  % logical array of true and false values
        jclassid = 'Z';
    case 'char'  % character array
        jclassid = 'C';
    case {'int8','uint8'}  % 8-bit signed and unsigned integer array
        jclassid = 'B';
    case {'int16','uint16'}  % 16-bit signed and unsigned integer array
        jclassid = 'S';
    case {'int32','uint32'}  % 32-bit signed and unsigned integer array
        jclassid = 'I';
    case {'int64','uint64'}  % 64-bit signed and unsigned integer array
        jclassid = 'J';
    case 'single'  % single-precision floating-point number array
        jclassid = 'F';
    case 'double'  % double-precision floating-point number array
        jclassid = 'D';
    case 'complexsingle'
        jclassid = 'Lhu.bme.aut.matlab.ComplexF;';
    case 'complexdouble'
        jclassid = 'Lhu.bme.aut.matlab.Complex;';
    otherwise
        error('java:javaclass:InvalidArgumentValue', ...
            'MatLab type "%s" is not supported in Java.', mtype);
end
jclassname = [repmat('[',1,ndims), jclassid];