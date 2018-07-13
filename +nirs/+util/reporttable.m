function cfr_table=reporttable(tbl)

cfr_table = rptgen.cfr_ext_table('IsPgwide',false,'NumCols',num2str(length(tbl.Properties.VariableNames)));

for idx=1:length(tbl.Properties.VariableNames)
    tblcolumn{idx}=rptgen.cfr_ext_table_colspec('ColNum',num2str(idx));
    tblcolumn{idx}.ColName=tbl.Properties.VariableNames{idx};
    setParent(tblcolumn{idx}, cfr_table);
end

cfr_table_head = rptgen.cfr_ext_table_head;
setParent(cfr_table_head,cfr_table);
cfr_head_row = rptgen.cfr_ext_table_row;
setParent(cfr_head_row,cfr_table_head);

for idx=1:length(tbl.Properties.VariableNames)
    
    cfr_head_entry(idx)=rptgen.cfr_ext_table_entry;
    cfr_head_paragraph(idx)=rptgen.cfr_paragraph;
    cfr_head_text(idx)=rptgen.cfr_text('Content',tbl.Properties.VariableNames{idx});
    
    set(cfr_head_paragraph(idx),'ParaTextComp', cfr_head_text(idx));
    setParent(cfr_head_paragraph(idx),cfr_head_entry(idx));
    setParent(cfr_head_entry(idx),cfr_head_row);
end



cfr_table_body = rptgen.cfr_ext_table_body;
setParent(cfr_table_body,cfr_table);
cnt=1;
for idx=1:height(tbl)
    cfr_table_row{idx} = rptgen.cfr_ext_table_row;
    setParent(cfr_table_row{idx},cfr_table_body);
    for idx2=1:length(tbl.Properties.VariableNames)
        cfr_ext_table_entry{cnt} = rptgen.cfr_ext_table_entry;
        cfr_paragraph{cnt} = rptgen.cfr_paragraph;
        
        entry=tbl.(tbl.Properties.VariableNames{idx2})(idx,:);
        if(isa(entry,'rptgen.cfr_image') | isa(entry,'rptgen.cfr_link'))
            setParent(entry,cfr_paragraph{cnt});
         
        else
            if(iscell(entry)); entry=entry{1}; end;
            if(isnumeric(entry)); entry=num2str(entry); end;
            
            
            cfr_text{cnt} = rptgen.cfr_text('Content',entry);
            set(cfr_paragraph{cnt},'ParaTextComp',cfr_text{cnt});
        end
        setParent(cfr_paragraph{cnt},cfr_ext_table_entry{cnt});
        setParent(cfr_ext_table_entry{cnt},cfr_table_row{idx});
        
        cnt=cnt+1;
    end
    setParent(cfr_table_row{idx},cfr_table_body);
    
end

