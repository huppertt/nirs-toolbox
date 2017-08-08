function algorithm = kmedoidsDetermineAlgorithm(X,k)
%kmedoidsDetermineAlgorithm - determine which algorithm to use in kmedoids
%partional clustering
%
% algorithm = kmedoidsDetermineAlgorithm(X,k)
% where,
% X is an N by P array of data supplied to kmedoids
% k is the number of clusters sought
%
% algorithm is a string indicating the algorithm to use.
%
%   The algorithm chosen by default depends on the size of X:
%       - If the number of rows of X is less than 3000, PAM is chosen. 
%       - If the number of rows is between 3000 and 10000, small is chosen. 
%       - For all other cases, large is chosen.
%   These defaults are chosen to give a reasonable balance between the time
%   to run the algorithm and the quality of the resulting solution. The
%   best choice will vary depending on the needs of the user and the data
%   being clustered.

% Copyright MathWorks 2014

if size(X,1)<3000
    algorithm = 'pam';
elseif size(X,1)<10000
    algorithm = 'small';
else
    algorithm = 'large';
end
    