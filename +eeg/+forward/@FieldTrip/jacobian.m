function J = jacobian( obj)

    cfg=obj.getFieldTripCFG;
    leadfield = ft_prepare_leadfield(cfg);
    L=[];
    for idx=1:length(leadfield.leadfield)
        n=size(leadfield.leadfield{idx},1);
        L(:,idx)=sqrt(sum(leadfield.leadfield{idx}.^2,2));
    end
    
    L(find(isnan(L)))=0;
    J.eeg=L;
    
end
