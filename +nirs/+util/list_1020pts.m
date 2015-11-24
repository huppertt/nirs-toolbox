function varargout = list_1020pts(str,headsize)

fid=fopen(which('ext1020.sfp'),'r');
marker=textscan(fid,'%s\t%d\t%d\t%d');
fclose(fid);

XYZ(:,1)=double(marker{2});
XYZ(:,2)=double(marker{3});
XYZ(:,3)=double(marker{4});

marker{1}(ismember(marker{1},'spmlpa'))=cellstr('lpa');
marker{1}(ismember(marker{1},'spmrpa'))=cellstr('rpa');
marker{1}(ismember(marker{1},'spmnas'))=cellstr('nas');

if(nargin==0 || any(ismember(str,'?')))
    lst=1:length(marker{1});
else
    lst=find(ismember(lower(marker{1}),lower(str)));
end

Name = marker{1}(lst);
X=XYZ(lst,1);
Y=XYZ(lst,2);
Z=XYZ(lst,3);

Type=repmat({'10-20'},size(X));
Units=repmat({'mm'},size(X));

tbl=table(Name,X,Y,Z,Type,Units);


%% Resize the 10-20 system to match the head size (if provided)
if(nargin>1)
    % First, find the default size of our head
    tbl = nirs.util.register_headsize(tbl,headsize);
end


if(nargout==0)
    disp(tbl);
else
    varargout{1}=tbl;
end