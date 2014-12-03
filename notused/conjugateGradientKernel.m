function [p, sse, covp, ci] = conjugateGradientKernel( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 30;
    end
    
    if nargin > 3
        tolX = varargin{2};
    else
        tolX = .001;
    end
    
    % initialize
    p = p0; d = zeros(size(p)); dp = Inf*p;
    
    [dy, J] = fun( p );
    sse = dy'*dy;
    gradF = -2*J'*dy;
    
    
    iter = 1;
    while iter == 1 || (any( abs(dp./p) > tolX ) && iter < maxIter)
        gradF0 = gradF;
        gradF = -2*J'*dy;
        
        % conjugate direction via fletcher-reeves
        b = (gradF'*gradF)/(gradF0'*gradF0);
        d = -gradF + b*d;
        
        dp = - sse/(gradF'*d) * d;
        dy = Inf;
        
        % bisection line search
        dp = dp*5;
        while dy'*dy > sse || all( abs(dp./p) < tolX )
            dp = dp/5;
            [dy, J] = fun( p + dp );
        end
        
        sse = dy'*dy;
        p = p + dp;   
        
        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
                
        iter = iter + 1;
    end
    
    v = (length(dy)-length(p));
    covp = pinv(J'*J) * sse / v ;
    delta = tinv( .975, v ) * sqrt( diag( covp ) );
    ci = [p-delta p+delta];
    
    
    
    
%     % initialize
%     p = p0;
%     dp = Inf;
%     iter = 1;
% 
%     [dy, J] = fun( p );
%         
%     % sum squared errors
%     F1 = dy'*dy;
% 
%     % gradient 
%     gradF1 = -2*J'*dy;
%     
%     
%     while any( abs(dp./p0) > tolX ) && iter < maxIter 
%         
%         disp( ['SSE: ' num2str(F1) ' Parameters: ' num2str( p' )] ); 
%         
%         % solve to find step size: 
%         % F(x+dx) = F(x) + gradF * (a*d) = 0
%         d = -gradF1;% ./ abs(gradF1) .* rand( size(gradF1) );
%         a = - F1 / (gradF1'*d);
%         
%         % check a
%         [dy, J] = fun( p + a*d );
%         
%         % update sse and gradient
%         F2 = dy'*dy;
%         gradF2 = -2*J'*dy;
%         
%         disp( ['SSE: ' num2str(F2) ' Parameters: ' num2str( (p + a*d)' )] ); 
%         
%         % reduction ratio
%         RR = (F1 - F2) / (-a*gradF1'*d); %% CHANGE TO GOLDSTEIN OR WOLFE
%         
%         if 1; %RR > .5 || gradF2'*d < gradF1'*d % concave down
%             % accept gradient descent step
%             p = p + a*d;
%             F1 = F2;
%             gradF1 = gradF2;
%         else
%             
% %         if gradF2*d < gradF1*d % concave down
% %             F1 = F2;
% %             gradF1 = gradF2;
% %                 
% %                 [dy, J] = fun( p + a*d );
% %                 
% %                 F2 = dy'*dy;
% %                 gradF2 = -2*J'*dy;
%                 
%             
%             % fit quadratic to two points
%             scal = abs(gradF1'*d);
%             a = a*scal;
%             
% %             A = [1 0 0 0;
% %                 0 1 0 0;
% %                 1 a a^2 a^3;
% %                 0 1 2*a 3*a^2];
% %             
% %             b = [F1 gradF1'*d/scal F2 gradF2'*d/scal]';
%             
%             A = [1 0 0;
%                 0 1 0;
%                 1 a a^2;
%                 0 1 2*a];
%             
%             b = [F1 gradF1'*d/scal F2 gradF2'*d/scal]';
% %             
%             c = double( pinv(A)*b );
%             
% %             f = @(a,c) c(1) + c(2)*a + c(3)*a.^2 + c(4)*a.^3;
% %             df = @(a,c) c(2)*a + 2*c(3)*a + 3*c(4)*a.^2;
% %             g = @(c) double( [(F1-f(0,c));
% %                 (gradF1'*d/scal-df(0,c));
% %                 (F2-f(a,c));
% %                 (gradF2'*d/scal-df(a,c))] );
% %             
% %             lsqnonlin( g, [0 0 0 0]',[-Inf -Inf -Inf 0]',[Inf Inf Inf Inf]' )
%             
% %             a = fsolve( @(a)c(2)+2*c(3)*a+3*c(4)*a.^2, double(a) );
%             a = -c(2)/c(3)/2;
%             
%             a = a/scal;
%             p = p + a*d;
%             
%             [dy, J] = fun( p );
%             F1 = dy'*dy;
%             gradF1 = -2*J'*dy;
%         end
%         
%         dp = a*d;
% 
%         % useful output
% %         disp( ['SSE: ' num2str(F1) ' Parameters: ' num2str( p' )] ); 
%         iter = iter + 1;
%     end
    
end