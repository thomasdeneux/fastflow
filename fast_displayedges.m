function [hlines, hnumbers] = fast_displayedges(CS,edges,cmd)
% function fast_displayedges(CS,edges,cmd)
%---
% displays edges with numbers and colors (edges should have a field
% 'points2' or a field 'snake')
% executes cmd when an edge is selected; cmd can be either a character
% array or a function with 2 parameters (i,ei)
% if cmd is a character array, use
% - i to design edge number
% - ei to design edges(i)

if nargin==0, help fast_displayedges, return, end

if nargin<3
    fun = '';
    ht = 'off';
else
    if ischar(cmd)
        str = cmd;
        cmd = @(i,ei)callback(i,ei,str);
    end
    fun = @(hl,evnt)cmd(getappdata(hl,'i'),getappdata(hl,'e'));
    ht = 'on';
end
    
colormap gray, imagesc(CS'), 
% if cbflag, fn_imvalue end, axis image, fn_imvalue, else axis image, end %#ok<DUALC>
nedges = length(edges);
for i=1:nedges
    if isfield(edges,'points2')
        x = edges(i).points2;
    elseif isfield(edges,'snake')
        x = edges(i).snake;
    else
        error('no field ''points2'' or ''snake''')
    end
    if size(x,1)>2, x=x'; end
    if isfield(edges,'color'), col = edges(i).color; else col = rand(1,3); end
    hlines(i) = line(x(1,:),x(2,:),'color',col,'linewidth',2,'hittest',ht,'buttondownfcn',fun);
    setappdata(hlines(i),'i',i), setappdata(hlines(i),'e',edges(i))
    hnumbers(i) = text('position',x([1 2],1),'color',col,'fontweight','bold','string',num2str(i));
end

if nargout==0, clear hlines hnumbers, end

%---
function callback(i,ei,str)

assignin('base','i',i)
assignin('base','ei',ei)
evalin('base',str)

