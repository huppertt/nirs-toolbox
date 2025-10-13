function varargout =toolbox_report_chapter()
% This function will report all the info for the version of matlab,
% nirs-toolbox, etc

warning('off','MATLAB:strrep:InvalidInputType');
matlabInfo=struct2table(ver);
info=nirs.get_toolbox_info;

matlabInfo(ismember(matlabInfo.Name,'+NIRS TOOLBOX'),:).Version{1}=info.Latest_Tag;
matlabInfo(ismember(matlabInfo.Name,'+NIRS TOOLBOX'),:).Date{1}=info.Commit_Date;
matlabInfo(ismember(matlabInfo.Name,'+NIRS TOOLBOX'),:).Release{1}=['https://www.github.com/huppertt/nirs-toolbox.git ' info.Commit_Hash_short];


import mlreportgen.report.*
import mlreportgen.dom.*

rpt_cpt=Chapter('Toolbox Version Info');

tbl=Table(matlabInfo);

tbl.Style = [tbl.Style
    {NumberFormat("%1.4f"),...
    Width("100%"),...
    Border("solid"),...
    ColSep("solid"),...
    RowSep("solid")}];
rpt_tbl=BaseTable(tbl);

rpt_cpt.add(rpt_tbl);

varargout{1}=rpt_cpt;

if(nargout>1)
    tbls_out.name='Toolbox information';
    tbls_out.table=matlabInfo;
    savedOutputs.tables(1)=tbls_out;
    varargout{2}=savedOutputs;
end