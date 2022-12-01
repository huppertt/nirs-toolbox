function DPF = frontal_DPF_model(lambda,data,Age)
%
% Felix Scholkmann, Martin Wolf,
% "General equation for the differential pathlength factor of the frontal human head depending on wavelength and age,"
% J. Biomed. Opt. 18(10) 105004 (11 October 2013) https://doi.org/10.1117/1.JBO.18.10.105004

persistent doonce;
if(nargin==3)
    if(isa(Age,'char') | isa(Age,'string'))
        if(data.demographics.iskey(Age))
            Age=data.demographics(Age);
        else
            if(isempty(doonce) || ~doonce)
                warning(['Demographics field: ' Age  ' not found']);
                doonce=true;
            end
            Age=21;
        end
    end
else
    Age=data;
end

if(nargin<2)
    Age=21;
end


alpha=223.3;
beta=0.05624;
gamma=0.8493;
delta= -5.723E-7;
epsilon= 0.001245;
epsi= -0.9025;

% Equation 7
DPF = alpha + beta*Age^gamma + delta*lambda.^3 + epsilon * lambda.^2 + epsi*lambda;