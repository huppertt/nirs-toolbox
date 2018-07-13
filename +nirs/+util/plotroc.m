function out1 = plotroc(varargin)
%PLOTROC Plot receiver operating characteristic.


hold on;
names={};
out1=[];
colors = lines(nargin/3);

cnt=1;
for i=1:3:nargin
     [tp,fp,th] = nirs.testing.roc(varargin{i},varargin{i+1});
     names{cnt}=varargin{i+2};
     out1(cnt)=plot(fp, tp, 'Color', colors(cnt,:));
     cnt=cnt+1;
end
legend(names, 'Location', 'SouthEast')
xlabel('False Positive Rate')
ylabel('True Positive Rate')