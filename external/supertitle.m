function supertitle(varargin)
if nargin<2
    fig=gcf; str=varargin{1};
else
    fig=varargin{1}; str=varargin{2};
end
str=strrep(strrep(str,'_','\_'),'\\_','\_');
allaxes=findall(fig,'type','axes');
allpos = zeros(length(allaxes),4);
for x=1:length(allaxes)
   allpos(x,:)=get(allaxes(x),'Position'); 
end

pause(.1);
shiftAxes(fig,'down',.05);
mtit(fig,sprintf('%s\n',str));

end
