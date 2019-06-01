function StartFastFlow(varargin)

% Add relevant folders to the path
basedir = fileparts(which('StartFastFlow'));
if exist(fullfile(basedir,'madflow.m'), 'file')
    % we are executing a 'StartFastFlow' that is inside the fastflow folder
    % -> basedir is one step above
    basedir = fileparts(basedir);
elseif exist(fullfile(basedir,'brick'), 'dir')
    % we are executing a 'StartFastFlow' that is above the brick, colormaps
    % and fastflow folders
    addpath(basedir)
else
    warning('Failed to locate brick, colormaps and fastflow folders')
end
addpath(fullfile(basedir,'brick'))
addpath(fullfile(basedir,'colormaps'))
addpath(fullfile(basedir,'fastflow'))

% Try to compile MEX files for movie coregistration if necessary
if ~exist('fast_regenergy','file')
    disp('COMPILE FILE fast_regenergy.cpp')
    swd = pwd;
    fil = which('fast_regenergy.cpp');
    p = fileparts(fil);
    fprintf('cd %s\n',p)
    cd(p)
    disp('mex fast_regenergy.cpp')
    try
        mex fast_regenergy.cpp
        fprintf('success\n')
    catch %#ok<*CTCH>
        fprintf('failed\n')
    end
    cd(swd)
end
if ~exist('fast_regenergy_single','file')
    disp('COMPILE FILE fast_regenergy_single.cpp')
    swd = pwd;
    fil = which('fast_regenergy_single.cpp');
    p = fileparts(fil);
    fprintf('cd %s\n',p)
    cd(p)
    disp('mex fast_regenergy_single.cpp')
    try
        mex fast_regenergy_single.cpp
        fprintf('success\n')
    catch
        fprintf('failed\n')
    end
    cd(swd)
end


% Start the program
madflow(varargin{:});

