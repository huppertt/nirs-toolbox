function out = parseTest()

f = 'y ~ 1 + a*b*c + (1|a)';

out = parseFormula(f);



end


function out = parseFormula( f )


    s = splitOn( f, '~' );
    
    lhs = s.left;
    rhs = s.right;
    
    s = strtrim(strsplit(rhs,'+'));
    
    out = splitOn( out, '+' );
    out = splitOn( out, '*' );

end


function out = splitOn( s, delim )

    if ischar(s)
        i = strfind(s,delim);

        if isempty(i)
            out = s;
        else
            i = i(1);

            out.op      = delim;
            out.left   	= strtrim( s(1:i-1) );
            out.right 	= splitOn( strtrim(s(i+1:end)), delim );
        end
    else
        out = s;
        out.left    = splitOn( out.left, delim );
        out.right   = splitOn( out.right, delim );
    end

end


%%
y = randn(100,1);
a = categorical(randi(3,[100 1]));
b = categorical(randi(3,[100 1]));
c = randn(100,1);
d = randn(100,1);

tbl = table(y,a,b,c,d);

lm = fitlme(tbl,'y ~ a*b*c*d + (1 | a)');

