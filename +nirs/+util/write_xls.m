function write_xls(filename,Data,sheetname)

root=fileparts(which('nirs.util.write_xls'));

%% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
if exist('org.apache.poi.ss.usermodel.WorkbookFactory', 'class') ~= 8 ...
        || exist('org.apache.poi.hssf.usermodel.HSSFWorkbook', 'class') ~= 8 ...
        || exist('org.apache.poi.xssf.usermodel.XSSFWorkbook', 'class') ~= 8
    javaaddpath(fullfile(root,'private','poi-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','poi-ooxml-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','poi-ooxml-schemas-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','xmlbeans-2.3.0.jar'));
    javaaddpath(fullfile(root,'private','dom4j-1.6.1.jar'));
    javaaddpath(fullfile(root,'private','stax-api-1.0.1.jar'));
end


if(isa(Data,'table'))
    names=Data.Properties.VariableNames;
    Data=table2cell(Data);
    Data=vertcat(names,Data);
end

if(nargin>2 && ~isempty(sheetname))
    warning('off','xlwrite:AddSheet');
    xlwrite(filename,Data,sheetname);
else
    xlwrite(filename,Data);
end