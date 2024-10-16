function [info,tbl]=load_BIDS_JSON(filename)

rootname=strtok(filename,'.');
jsonfile=[rootname '.json'];
tsvfile=[rootname '.tsv'];

info=textread(jsonfile,'%s');

if(info{1}=='{' & info{2}(1)=='{')
    info={info{2:end-1}};
end
info=strcat(info{:});
info=strrep(info,'""','""');
info=strrep(info,'},{',',');

info=jsondecode(info);

if(exist(tsvfile))
    tbl=readtable(tsvfile,'FileType','text');
else
    tbl=[];
end