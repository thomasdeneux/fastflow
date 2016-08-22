function fastflow_package

%% Version

version = '2.2.1'; 


%% Create empty base folder

disp 'create empty base folder'
matlabdir = fn_cd('matlab');
basedir = fn_cd('matlab',['FastFlow' '-' version]);
if exist(basedir,'dir')
    try
        rmdir(basedir,'s')
    catch ME
        warning(ME.message)
    end
end
mkdir(basedir)

%% Encrypt m-files and copy regular files

% StartFastFlow
disp 'StartFastFlow.m'
cd(fullfile(matlabdir,'fastflow'))
copyfile('StartFastFlow.m',basedir)

% brick
folder = 'brick';
mfiles = listfiles(folder,'*.m',{'fn_dodebug.m' 'fn_cd.m'});
deployregular(mfiles,folder,basedir)
disp(fullfile(folder,'fn_dodebug.m'))
fn_savetext({'function b = fn_dodebug','','b = false;'},'fn_dodebug.m')
regularfiles = listfiles(folder,'*.png');
deployregular(regularfiles,folder,basedir)

folder = fullfile('brick','private');
mfiles = listfiles(folder,'*.m');
deployregular(mfiles,folder,basedir)

% colormaps
folder = 'colormaps';
mfiles = listfiles(folder,'*.m');
deployregular(mfiles,folder,basedir)

% fastflow
folder = 'fastflow';
mfiles = listfiles(folder,'*.m',{'StartFastFlow.m' 'fastflow_package.m'});
deployregular(mfiles,folder,basedir)
regularfiles = {'madflow.mat' 'fast_color.mat' 'fast_regenergy.cpp' 'fast_regenergy_single.cpp'};
deployregular(regularfiles,folder,basedir)


%% License and doc

folder = fullfile('fastflow','license and doc');
regularfiles = {'license.pdf' 'license.txt' 'README.txt'};
deployregular(regularfiles,folder,basedir,'')

%% Version
fn_savetext({['This is fastflow version ' version]}, ...
    fullfile(basedir,['fastflow version ' version '.txt']))

%---
function b = copyhelp(fun,target)
% function copyhelp(fun,target)
%---
% copy only the help content of an m-file

% get the help content
fun = fn_fileparts(fun,'base');
txt = help(fun);
if isempty(txt), b = false; return; end
txtLineByLine = strcat('%',regexp(txt,'\n','split')');

% build the target file name
if exist(target,'dir')
    filename = fullfile(target,[fun '.m']);
else
    filename = target;
end

% create the m-file with help content
fid = fopen(filename,'w');
if (fid == -1)
    disp('Unable to write to %s', filename);
end
fprintf(fid, '%s\n',txtLineByLine{:});
fclose(fid);
b = true;

%---
function files = listfiles(folder,pattern,excludelist)

source = fn_cd('matlab',folder);
files = dir(fullfile(source,pattern));
files = {files.name};
if nargin>=3
    files = setdiff(files,excludelist);
end

%---
function deploypcode(mfiles,folder,basedir,targetfolder)

source = fn_cd('matlab',folder);
if nargin<4, targetfolder = folder; end
target = fullfile(basedir,targetfolder);

mkdir(target)
cd(target)
for k=1:length(mfiles)
    fk = mfiles{k};
    if copyhelp(fk,'.')
        disp(fullfile(targetfolder,fk))    
    end
    disp(fullfile(targetfolder,[fk(1:end-2) '.p']))    
    pcode(fullfile(source,fk))
end

%---
function deployregular(regularfiles,folder,basedir,targetfolder)

source = fn_cd('matlab',folder);
if nargin<4, targetfolder = folder; end
target = fullfile(basedir,targetfolder);

if ~exist(target,'dir'), mkdir(target), end
cd(target)
for k=1:length(regularfiles)
    fk = regularfiles{k};
    disp(fullfile(targetfolder,fk))    
    copyfile(fullfile(source,fk),'.')
end