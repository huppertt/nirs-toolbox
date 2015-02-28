function flag = canwrite(loc)

% flag = canwrite(loc)
%
% determines if there is write access to the directory
%
% loc is a location on the computer. If no loc is given,
%   the current Matlab directory will be checked
% flag is the result (1 if there is write access, 0 otherwise)
% alt_path is the path that user have write access to.
% 

if nargin == 0 || isempty(loc)
    loc = pwd;
end

if isdir(loc)
    fooloc = loc;
else
    if isempty(fileparts(loc))
        fooloc = pwd;
    else
        fooloc = fileparts(loc);
    end
end

% k1 = strfind(loc,filesep);
% 
% if isempty(k1)
%     loc = pwd;
% end

% if ~isdir(loc)
%     loc = fileparts(loc);
% end

[foo1 testFile ext] = fileparts(tempname);
loc = fooloc;

if exist(loc, 'dir') == 7
    [fid message] = fopen([loc filesep testFile ext],'wt');
    if fid ~= -1 && isempty(message)
        fclose(fid);
        delete([loc filesep testFile ext]);
        flag = 1;
    else
        flag = 0;
    end
else
    % loc doesn't exist
    flag = 0;
end