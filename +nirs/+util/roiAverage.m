function tbl = roiAverage( data, R, names )
    if ischar( names )
        names = {names};
    end
    
     % sort probe
    link = data.probe.link;
    [link, ilink] = sortrows(link, {'type','source', 'detector'});
    
    if(isa(R{1},'table'))
        
        %First deal with the NaN values;
        allSrc=unique(link.source);
        allDet=unique(link.detector);
        
        for idx=1:length(R)
            src=R{idx}.source(isnan(R{idx}.detector));
            if(length(src)>0)
            R{idx}=[R{idx}; table(kron(src,ones(length(allDet),1)),...
                kron(ones(length(src),1),allDet),'VariableNames',{'source','detector'})];
            end
            
            det=R{idx}.detector(isnan(R{idx}.source));
            if(length(det)>0)
            R{idx}=[R{idx}; table(kron(ones(length(det),1),allSrc),...
                    kron(ones(length(allSrc),1),det),'VariableNames',{'source','detector'})];
            end
            
        end
        
        
        % The region definition is a table, parse it to the contrast vector
        types=unique(link.type);
        
        RNew=cell(length(R)*length(types),1);
        NamesNew=cell(length(R)*length(types),1);
        cnt=1;
        for idx=1:length(types)
            for idx2=1:length(R)
                RNew{cnt}=ismember(link,[R{idx2} table(repmat({types{idx}},height(R{idx2}),1),...
                    'VariableNames',{'type'})]);
                NamesNew{cnt}=[names{idx2} ':' types{idx}];
                cnt=cnt+1;
            end
        end
        R=RNew;
        names=NamesNew;
    end
    
   
    
    % change ROIs to sorted indices
    for i = 1:length(R)
        R{i} = ilink(R{i});
    end
    
    % sort variables
    [vars, ivars] = sortrows(data.variables, {'cond', 'source', 'detector', 'type'});
    beta = data.beta(ivars);
    covb = data.covb(ivars, ivars);
    
    % unique conditions
    uconds = unique(vars.cond, 'stable');
    
    % loop over conditions
    varnames = {'ROI', 'Contrast', 'Beta', 'SE', 'DF', 'T', 'p'};
    tbl = table({},[],[],[],[],[],[], 'VariableNames', varnames);

    for i = 1:length(uconds)
        lst = strcmp(vars.cond, uconds{i});
        b = beta(lst);
        C = covb(lst,lst);
        
        for j = 1:length(R)
            % contrast vector
            c = zeros(size(b));
            c(R{j}) = 1/length(R{j});
            
            broi    = c'*b;
            se      = sqrt(c'*C*c);
            t       = broi / se;
            df      = data.dfe;
            p       = 2*tcdf(-abs(t),df);
            
            tmp = cell2table({names{j}, uconds{i}, broi, se, df, t, p});
            tmp.Properties.VariableNames = varnames;
            tbl = [tbl; tmp];
        end
    end
    
    q   = nirs.math.fdr( tbl.p );
    tbl = [tbl table(q)];
end