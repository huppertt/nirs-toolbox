function QCReport(data)

figure;
n=length(data);

for i=1:n;
    subplot(4,n,i);
    data(i).draw;
    legend off;
    
    sni=nirs.math.structnoiseindex(data(i).data);
    types=unique(data(i).probe.link.type);
    lst=ismember(data(i).probe.link.type,types(1));
    [colors,lineStyles]=nirs.util.values2lines(sni(lst),[-10 10]);
    data(i).probe.draw(colors,lineStyles,subplot(4,n,n+i));
      lst=ismember(data(i).probe.link.type,types(2));
    [colors,lineStyles]=nirs.util.values2lines(sni(lst),[-10 10]);
    data(i).probe.draw(colors,lineStyles,subplot(4,n,2*n+i));
    
    subplot(4,n,3*n+i);
    hist(sni);
end