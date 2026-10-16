function set_up_dirs(all)
%SET_UP_DIRS  Add required toolbox directories to the MATLAB path.
%
%   nirs.set_up_dirs()       — adds commonly needed external libraries
%   nirs.set_up_dirs(true)   — adds ALL external subdirectories via genpath

if nargin < 1
    all = false;
end

% Resolve toolbox root: this file lives in +nirs/, so go up one level
thisFile = mfilename('fullpath');
[pkgDir, ~, ~] = fileparts(thisFile);   % +nirs/
[rootDir, ~, ~] = fileparts(pkgDir);    % nirs-toolbox/

extDir = fullfile(rootDir, 'external');

if all
    addpath(genpath(extDir));
    addpath(genpath(fullfile(rootDir, 'demos')));
    fprintf('Added all external and demo paths.\n');
    return
end

% Selective: commonly needed externals
libs = { ...
    ''                      % external/ root (has standalone .m files)
    'Dictionary'
    'cbrewer'
    'nirsviewer'
    'export_fig'
    'aal'
    'fieldtrip'
    'spm'
    'nirfast'
    'view_nii'
    'fastica'
    'icp'
    'treeTable'
    'bvaloader'
    'edfreadZip'
    };

for i = 1:numel(libs)
    p = fullfile(extDir, libs{i});
    if isfolder(p)
        addpath(p);
    end
end

% Demos (small, add recursively)
demosDir = fullfile(rootDir, 'demos');
if isfolder(demosDir)
    addpath(genpath(demosDir));
end

fprintf('Paths added.\n');

end
