function [mdc,pwr] = MDC(Stats,power,alpha)

side=1;

beta=1-power;

if(isa(Stats,'nirs.core.ChannelStats'))    
    tbeta = tinv(1-beta,Stats.dfe-2);
    talpha = tinv(1-alpha/side,Stats.dfe-2);
    mdc = sqrt(diag(Stats.covb)).*(tbeta+talpha);
    
    %Stats.beta./(sqrt(diag(Stats.covb))=t2beta+talpha;
   
    talpha2 = tinv(1-alpha/side,Stats.dfe-2);
    pwr=1-tcdf((abs(Stats.tstat)-talpha2),Stats.dfe-2);
elseif(isa(Stats,'table'))
    t2beta = tinv(1-beta,Stats.DF-2);
    talpha = tinv(1-alpha/side,Stats.DF-2);
    mdc = Stats.SE.*(t2beta+talpha);
    talpha2 = tinv(1-alpha/side,Stats.DF-2);
    pwr=1-tcdf(abs(Stats.T)-talpha2,Stats.DF-2);
elseif(isa(Stats,'nirs.core.ImageStats'))
    t2beta = tinv(1-beta,Stats.dfe-2);
    talpha = tinv(1-alpha/side,Stats.dfe-2);
    
    c=diag(Stats.covb_chol*Stats.covb_chol')+Stats.typeII_StdE.^2;
    mdc = sqrt(c).*(t2beta+talpha);
    
    %Stats.beta./(sqrt(diag(Stats.covb))=t2beta+talpha;
    
    talpha2 = tinv(1-alpha/side,Stats.dfe-2);
    pwr=1-tcdf(abs(Stats.tstat)-talpha2,Stats.dfe-2);
end
