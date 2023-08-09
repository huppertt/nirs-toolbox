function [ind_big, int_func_big]=mytsearchn_bem(mesh,coord)

% [ind_big, int_func_big]=mytsearchn_bem(mesh,coord)
%
% calculates integrating functions for measurement position; 
% calculated by finding nearest node and element containing the measurement
% location; and using areas to compute the integrating function
%
% mesh is the input mesh
% coord is the set of points
% ind_big are the indices for the elements containing the points
% int_func_big are the barycentric coordinates

int_func_big = [];
ind_big = [];
meas=coord;
nodes=mesh.nodes;
elements=mesh.elements;

for i = 1:size(meas,1)
    
    [ind_nod] = dsearchn(nodes,meas(i,:));
    elem1 = find(elements(:,1) == ind_nod);
    elem2 = find(elements(:,2) == ind_nod);
    elem3 = find(elements(:,3) == ind_nod);
    ind_elem = unique([elem1; elem2; elem3]);
     count = 1;
    for ii = 1:length(ind_elem)
       
        nod1 = nodes(elements(ind_elem(ii),1),:);
        nod2 = nodes(elements(ind_elem(ii),2),:);
        nod3 = nodes(elements(ind_elem(ii),3),:);
       [area] = get_area(nod1,nod2,nod3);
       %%first basis function
       nod1 = meas(i,:);
       nod2 = nodes(elements(ind_elem(ii),2),:);
       nod3 = nodes(elements(ind_elem(ii),3),:);
       [area1] = get_area(nod1,nod2,nod3);
       %%second basis function
       nod1 = meas(i,:);
       nod2 = nodes(elements(ind_elem(ii),3),:);
       nod3 = nodes(elements(ind_elem(ii),1),:);
       [area2] = get_area(nod1,nod2,nod3);
       %%third basis function
       nod1 = meas(i,:);
       nod2 = nodes(elements(ind_elem(ii),1),:);
       nod3 = nodes(elements(ind_elem(ii),2),:);
       [area3] = get_area(nod1,nod2,nod3);
       int_func(1) = area1/area;
       int_func(2) = area2/area;
       int_func(3) = area3/area;
       if (sum(int_func) < 1.03)
           break;
       end
       temp(count,1:3) = int_func;
       jj = sum(temp,2);
       cc = find(jj == min(jj));
       if length(cc) > 1
           cc = cc(1);
       end
       int_func(1:3) = temp(cc,1:3);
       clear jj cc;
       count = count+1;
    end
    clear temp;
    %sum(int_func)
       %%to make sure integral is 1.0
       int_func(1:end) = 1./sum(int_func).*int_func(1:end);
     %disp('sum(int_func) = ');
     %disp(sum(int_func))
    int_func_big = [int_func_big; int_func];
    ind_big = [ind_big;ind_elem(ii)];
    clear int_func;
    clear ind_elem;
end



 function [delG] = get_area(nod1,nod2,nod3)
     
 %%assuming linear triangular element, delG gives the conversion going from
 %%global to local co-ordinates; 
 %%normal components for the normal to the element, are also calculated.
 
 delG = sqrt((((nod1(2) - nod3(2))*(nod2(3) - nod3(3))) - ((nod2(2) - nod3(2))*(nod1(3) - nod3(3))))^2 + ...
     (((nod1(3) - nod3(3))*(nod2(1) - nod3(1))) - ((nod1(1) - nod3(1))*(nod2(3) - nod3(3))))^2 + ...
     (((nod1(1) - nod3(1))*(nod2(2) - nod3(2))) - ((nod1(2) - nod3(2))*(nod2(1) - nod3(1))))^2);
 