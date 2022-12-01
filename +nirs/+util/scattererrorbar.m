function h=scattererrorbar(x,y,e)

hold on;

h=[];
a=(max(x)-min(x))/10;
if(a<.1); a=.1; end;
    
for i=1:length(x)
    h(end+1)=line([x(i) x(i)],[y(i)-e(i) y(i)+e(i)]);
    h(end+1)=line([x(i)-a x(i)+a],[y(i)-e(i) y(i)-e(i)]);
    h(end+1)=line([x(i)-a x(i)+a],[y(i)+e(i) y(i)+e(i)]);
end


return
