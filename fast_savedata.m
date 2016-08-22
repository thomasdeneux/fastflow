function fast_savedata(Y,fname,model)
%---
% function fast_savedata(Y,fname,model)

if nargin<1, help fast_savedata, return, end

disp('change working directory disabled - assumed datype = 0')
datatype = 0;
swd = pwd;
catflag = false;
% % change working directory
% swd = pwd;
% datatype = 0;
% switch fname(1:2)
%     case 'TC'
%         catflag = false;
%         datatype = 0;
%         switch fname(20)
%             case '1'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%             case '2'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%             case '3'
%                 cd /storm5/thomas/Data/0503_Marseille/TOM03oct05
%         end
%     case 'fl'
%         catflag = true;
%         datatype = 2;
%         cd C:\Users\Thomas\WORK\IRMf\0503_Marseille\data\catflow\
% end

% get headers
[nfrms,xs,ys,hsize,nstims]=get_head(model,'nframesperstim,framewidth,frameheight,lenheader,nstimuli');

% Input
% frames
frames = 1:nfrms;
nt = length(frames);
% sub-window
subx = 1:xs; suby = 1:ys;
indices = 1:xs*ys;

% save data
if exist('subx','var')
    y = reshape(Y,length(subx)*length(suby),nt);
end
copyfile(model,fname)
savesum2(y,fname,frames-1,0,datatype,xs,ys,hsize,indices);

% back to old working directory
cd(swd)