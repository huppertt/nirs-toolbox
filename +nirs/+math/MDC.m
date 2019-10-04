function [mdc,pwr] = MDC(Stats,beta,alpha)

side=1;

if(isa(Stats,'nirs.core.ChannelStats'))
    
    t2beta = tinv(beta,Stats.dfe-2);
    talpha = tinv(1-alpha/side,Stats.dfe-2);
    
    mdc = sqrt(diag(Stats.covb)).*(t2beta+talpha);
    
    %Stats.beta./(sqrt(diag(Stats.covb))=t2beta+talpha;
    
    talpha2 = tinv(1-max(Stats.p,alpha)/side,Stats.dfe-2);
    pwr=tcdf(abs(Stats.tstat)-talpha2-t2beta,Stats.dfe-2);
elseif(isa(Stats,'table'))
    t2beta = tinv(beta,Stats.DF-2);
    talpha = tinv(1-alpha/side,Stats.DF-2);
    mdc = Stats.SE.*(t2beta+talpha);
    talpha2 = tinv(1-max(Stats.p,alpha)/side,Stats.DF-2);
    pwr=tcdf(abs(Stats.T)-talpha2-t2beta,Stats.DF-2);
elseif(isa(Stats,'nirs.core.ImageStats'))
    t2beta = tinv(beta,Stats.dfe-2);
    talpha = tinv(1-alpha/side,Stats.dfe-2);
    
    c=diag(Stats.covb_chol*Stats.covb_chol')+Stats.typeII_StdE.^2;
    mdc = sqrt(c).*(t2beta+talpha);
    
    %Stats.beta./(sqrt(diag(Stats.covb))=t2beta+talpha;
    
    talpha2 = tinv(1-max(Stats.p,alpha)/side,Stats.dfe-2);
    pwr=tcdf(abs(Stats.tstat)-talpha2-t2beta,Stats.dfe-2);
end
