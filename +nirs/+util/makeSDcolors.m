function colors_out=makeSDcolors(link,cmap_function)
% Assign same color to each unique SD combo
% Example Optode 12 should have the color 'Red' regardless of which
% wavelength is being plotted

% Takes in link of SD pairs and returns color array of identical side with
% each type for an SD pair assigned the same color


link_no_type=link;
link_no_type.type=[];
[ulink,~,idx_b]=unique(link_no_type(:,[2,1]),'rows');
num_u_colors=height(ulink);

if(nargin<2)
    cmap_function='lines';
end

try
    cmap=eval(cmap_function,num_u_colors); %ex 'lines'
catch
    cmap=lines(num_u_colors);
end

colors_out=zeros(height(link),3);
%[~,type_color_idx]=ismember(link.type,type);

for id=1:num_u_colors
    assign_idx=(id==idx_b);
    colors_out(assign_idx,:)=repmat(cmap(id,:),[sum(assign_idx),1]);
end