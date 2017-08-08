function p = bvncdf(b,rho,~)
% CDF for the bivariate normal.
 
% Implements Section 2.4 of Genz (2004).

n = size(b,1);
if rho == 0
    p = cast(prod(Phi(b),2), superiorfloat(b,rho));
else
    if abs(rho) < 0.3      % 6 point Gauss Legendre abscissas and weights
        w = [0.4679139345726904  0.3607615730481384  0.1713244923791705];
        y = [0.2386191860831970  0.6612093864662647  0.9324695142031522];
    elseif abs(rho) < 0.75 % 12 point Gauss Legendre abscissas and weights
        w = [0.2491470458134029  0.2334925365383547  0.2031674267230659 ...
             0.1600783285433464  0.1069393259953183  0.04717533638651177];
        y = [0.1252334085114692  0.3678314989981802  0.5873179542866171 ...
             0.7699026741943050  0.9041172563704750  0.9815606342467191];
    else                 % 20 point Gauss Legendre abscissas and weights
        w = [0.1527533871307259  0.1491729864726037  0.1420961093183821  0.1316886384491766  0.1181945319615184 ...
             0.1019301198172404  0.08327674157670475 0.06267204833410906 0.04060142980038694 0.01761400713915212];
        y = [0.07652652113349733 0.2277858511416451  0.3737060887154196  0.5108670019508271  0.6360536807265150 ...
             0.7463319064601508  0.8391169718222188  0.9122344282513259  0.9639719272779138  0.9931285991850949];
    end
    
    if abs(rho) < .925
        p1 = prod(Phi(b),2);
        asinrho = asin(rho);
        w = fliplr([w w]);
        theta = asinrho .* (fliplr([y -y]) + 1) / 2;
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        h = -b(:,1); k = -b(:,2);
        hk = h .* k;
        ssq = (h.^2 + k.^2) / 2;
        p2 = zeros(size(p1),class(rho));
        for i = 1:n
            if isfinite(hk(i))
                f = exp( -(ssq(i) - hk(i)*sintheta) ./ cossqtheta );
                p2(i) = asinrho * sum(w.*f) / 2;
            else
                % This piece is zero if either limit is +/- infinity.  If
                % either is NaN, p1 will already be NaN.
            end
        end
        p = p1 + p2/(2*pi);

    else % abs(rho) >= .925
        if rho > 0
            p1 = Phi(min(b,[],2));
            p1(any(isnan(b),2)) = NaN;
        else
            p1 = Phi(b(:,1)) - Phi(-b(:,2));
            p1(p1<0) = 0; % max would drop NaNs
        end
        
        s = sign(rho);
        if abs(rho) < 1
            h = -b(:,1); k = -b(:,2);
            shk = s.*h.*k;
            asq = 1 - rho.^2;  a = sqrt(asq);
            b = abs(h - s.*k); bsq = b.^2;
            c = (4 - shk)/8;
            d = c .* (12 - shk)/16;

            t1 = a.*(1 + d.*asq.^2/5 + (c/3 - d.*bsq/15).*(asq - bsq));
            t2 = b.*(1 - c.*bsq/3 + d.*bsq.^2/15);
            p2 = exp(-shk/2) ...
                .* (t1.*exp(-bsq./(2*asq)) - t2.*sqrt(2*pi).*Phi(-b./a));

            w = [w w];
            x = a .* ([y -y] + 1) / 2; x2 = x.^2; x4 = x2.^2;
            sqrt1mx2 = sqrt(1 - x2);
            t = (sqrt1mx2 - 1) ./ (2*(sqrt1mx2 + 1));
            for i = 1:n
                if isfinite(shk(i))
                    f = exp(-(bsq(i)./x2 + shk(i))/2) ...
                        .* (exp(shk(i).*t)./sqrt1mx2 - (1 + c(i).*x2 + d(i).*x4));
                    p2(i) = p2(i) + a .* sum(w.*f) / 2;
                else
                    % This piece is zero if either limit is +/- infinity.  If
                    % either is NaN, p1 will already be NaN.
                    p2(i) = 0;
                end
            end
        else % abs(rho) == 1
            p2 = zeros(class(rho));
        end
        p = p1 - s.*p2/(2*pi);
    end
end

end % bvncdf

function p = Phi(z)
% CDF for the normal distribution.
p = 0.5 * erfc(-z / sqrt(2));
end