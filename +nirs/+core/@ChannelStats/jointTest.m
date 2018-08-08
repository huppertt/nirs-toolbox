function out = jointTest( obj )
    %% jointTest - performs a joint hypothesis test across all channels in
    %              each SD pair (i.e. joint test of hbo & hbr for S1-D1)
    
    % ensure that data is sorted
    obj = obj.sorted();
    
    % get unique sd pair/conditions to jointly test
    vars = obj.variables;
    vars.type = [];
    
    if(ismember('model',vars.Properties.VariableNames))
        vars.model=[];
    end
    
    [~, uvars ,utests] = unique(vars, 'rows', 'stable');
    
    % loop through pairs/conditions
    for i = 1:max(utests)
       m = utests == i;
       covb=obj.covb(m,m);
       
       % sometimes not positive/definant due to numerical issues. 
       [U,s,V]=svd(covb);
       covb=0.5*(U*s*U'+V'*s*V);
       
       T2 = obj.beta(m)'*pinv(covb)*obj.beta(m);
       
       n = obj.dfe;
       k = sum(m);
       
       %Saitterwaitte estimate of effective DOF
       n=trace(covb)^2/(trace(covb^2))/(1/(n-1));
       
       
       F(i,1) = (n-k+1) ./ k ./ n .* T2;
       df1(i,1) = k;
       df2(i,1) = n-k+1;
    end
    
    % output
    out = nirs.core.ChannelFStats();
    
    out.F   = F;
    out.df1 = df1;
    out.df2 = df2;
    
    % new variables
    newVars = obj.variables(uvars, :);
    newVars.type = cell(size(newVars.type));
    newVars.type(:) = {'joint'};
    
    out.variables = newVars;
    
    % new probe
    link = obj.probe.link;
    link.type = [];
    
    [~,idx] = unique(link, 'rows', 'stable');
    
    out.probe = obj.probe;
    out.probe.link = out.probe.link(idx, :);
    out.probe.link.type = cell(size(out.probe.link.type));
    out.probe.link.type(:) = {'joint'};
    S.description = 'Joint F-test';
end