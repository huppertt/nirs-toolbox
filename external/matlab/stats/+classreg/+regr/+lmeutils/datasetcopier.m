function ds = datasetcopier(ds,D)

%   Copyright 2012-2013 The MathWorks, Inc.    

% Copy stuff from D into ds
props = D.Properties;

vn = props.VarNames;
for j=1:length(vn)
    vnj = vn{j};
    ds.(vnj) = D.(vnj);
end

fn = fields(props);
for j=1:length(fn)
    fnj = fn{j};
    ds.Properties.(fnj) = D.Properties.(fnj);
end
