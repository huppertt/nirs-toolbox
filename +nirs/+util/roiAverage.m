function tbl = roiAverage( data, R, names )
        if ischar( names )
            names = {names};
        end

        varnames = {'ROI', 'Contrast', 'Beta', 'SE', 'DF', 'T', 'p', 'q'};
        tbl = table({},[],[],[],[],[],[], 'VariableNames', varnames);

        R = R > 0;

        for i = 1:size(R,1)
            for j = 1:size(data.beta, 1)
                roi     = names{i};
                con     = data.names{j};

                w       = 1./squeeze( data.covb(j,j,R(i,:)) );
                x       = ones(size(w));

                ix      = pinv( x'*diag(w)*x ) * x' * diag(w);

                beta    = ix * data.beta(j,R(i,:))';
                se      = sqrt( pinv( x'*diag(w)*x ) );


                df      = data.dfe;
                t       = beta / se;
                p       = 2*tcdf(-abs(t),df);
                
                q       = nirs.math.fdr( p );

                tmp = cell2table({roi, con, beta, se, df, t, p, q});
                tmp.Properties.VariableNames = varnames;
                tbl = [tbl; tmp];
            end
        end
    end