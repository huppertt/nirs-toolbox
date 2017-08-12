function data=get_boundary_data_bem(elements,meas_int_func,phi)
% Calculates the phase and amplitude and location of detectors specified by
% 'meas_int_func' by using the solution values from 'phi'

nrow = size(meas_int_func,1);
data=zeros(1,nrow);
for j = 1 : nrow
    vtx_ind = elements(meas_int_func(j,1),:);
    data(j) = meas_int_func(j,2:end)*phi(vtx_ind);
end
