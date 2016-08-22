function StartFastFlow(varargin)

% Add relevant folders to the path
basedir = fileparts(which('StartFastFlow'));
addpath(basedir)
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

