function x = fast_loadfastflow(fname)
% function [S =] fast_loadfastflow(fname)

load(fname) % creates fastflow object 'x'

swd = pwd; cd(fileparts(fname)), rep = [pwd '/']; cd(swd)
x.savedir = rep;

if nargout==0
    assignin('base','S',x)
    clear x
end
