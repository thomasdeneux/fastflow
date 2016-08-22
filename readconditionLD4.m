% [timeaxis,conditions,multiplicity,trnum] = readconditionLD4(infile,stimlist,channels);
% 
% reads 7th data channel from REX B-files. 
% CAREFUL: THE STIMLIST MUST CONTAIN ALL THE STIMULI AND MUST BE CORRECT; OTHERWISE THE RESULT
% IS INCORRECT;
% Input: 
%   INFILE: name of file to read (string, e.g. 'myfileB'). Must be a B-type file. 
%   STIMLIST: Ecodes of conditions to read. Normally between 4000 and 8000. Row vector.
%   CHANNELS: Which channels to keep in conditions [laser doppler is number
%   7 usually]
%   
% Output:
%   TIMEAXIS: timeaxis, in ms, centered on stimulus onset. Row vector.
%   CONDITIONS: 5th analog REX-channel for the specified stimuli. 1st dim: data points (time), 
%   2nd dim: 5th  column of the B file: "doppler" (the different trials, without the rejected ones). Sequence: first 
%   all trials of one condition, then all of the next etc. Condition-ID order as in STIMLIST.
%   Data are cut to the maximum time interval common to all trials, between opening and closing daq window, 
%   centered on stimbegin.
%   MULTIPLICITY: Number of trials for each stimulus. Same ordering as in
%   stimlist. Does not assume equal number of trials per condition.
%   TRNUM: for each trial in CONDITIONS, its number in the acquisition
%   order (rejected trials not included)
% 
% IMPORTANT: Assumes that the code for the rejection is given twice. Also
% assumes that the file has LDind columns in standard order (time, ecodes, 4eyedata, doppler)+5 dummiy columns.
% The following typical sequence (with a rejection) shows the assumptions
% about the employed Ecode conventions and the organization of a trial
% 1001 trial begin
% stimID ecode of stimulus
% 800 open daq window
% 802 reject trial
% 802 reject trial (bis)
% 1001 trial begin
% stimID ecode of stimulus
% 800 open daq window
% 3232 begin stim
% 802 reject trial
% 802 reject trial (bis)
% 1001 begin trial 
% 800 open daq window
% 3232 begin stim
% 801 close daq window
% It is also assumed that there are no lost ms between daqwin opening and
% closing. 
%
% see also: READCONDITIONCAL, EXTRACTSTIMS, convertBfiles,
% splineminsquare, killoutliers2d, eyeplot
%
% Ivo, 1 july 2003, last revised 9 july 2003

function [timeaxis,conditions,multiplicity,trnum]=readconditionLD4(infile,stimlist,LDind);

if nargin<3
    LDind = 7;
end
LDind = LDind-2;
if any(LDind<3) || any(LDind>5) || any(mod(LDind,1)), error('condition channels have to be among 5,6,7'), end

% READING IN THE DATA
disp('reading in the data')
if ischar(infile)
    s = readconditionfile(infile);
else
    s = infile;
end
% check for LDind columns only
if size(s,2)~=5
    error('data should have 5 columns : time - flags - oxygen or junk - digital square - laser doppler (columns 1,2,5,6,7 from REX file)')
end

% LOOKING FOR INDICES
disp('looking for indices')
% setting the time sequences for each trial. CAREFUL! THESE ARE column INDICES of s, not times.

trialbegind=find(s(:,2)==1001); pause(0)
stimIDind=find(ismember(s(:,2),stimlist)); pause(0)
winopenind=find(s(:,2)==800); pause(0)
stimbegind=find(s(:,2)==3232); pause(0)
stimendind=find(s(:,2)==3333); pause(0)
rejectind=find(s(:,2)==802); pause(0) % they appear twice, one after the other.
wincloseind=find(s(:,2)==801); pause(0)

%creating errormessages if the files are "irregular"
if (max(size(trialbegind))~=max(size(winopenind)))
    error('not same number of trialbegins and data window openings in file')
elseif (max(size(trialbegind))~=(floor(max(size(rejectind))/2)+max(size(wincloseind))))
    error('not same number of trialbegins and (rejections+data window closings) in file')
end

% typical sequence with rejection:
% 1001 trial begin
% stimID ecode of stimulus
% 800 open daq window
% 802 reject trial
% 802 reject trial (bis; maybe to close daq window ???)
% 1001 trial begin
% stimID ecode of stimulus
% 800 open daq window
% 3232 begin stim
% 802 reject trial
% 802 reject trial 
% 1001 begin trial 
% 800 open daq window
% 3232 begin stim
% 801 close daq window

%-------------------------
% REMOVING REJECTED TRIALS
%-------------------------
disp('removing rejected trials')
if (isempty(rejectind) ~=1)
   % deletion of the - eventually existing - rejected trals

   % finding the events (trialbegins and data-window openings) corresponding to rejected trials. The if loop is to
   % take care of situations where the ecode nearest to a rejection-ecode is
   % belonging to the next good trial instead than to the one to be rejected.
   % these 2 vectors now contain the INDICES corresponding to the rejected trials (not their elements!!).
   indtrialbegdel=[];
   indwinopendel=[];
   for i = 1:2:max(size(rejectind)) 
      [m dummy] = max(find(s(rejectind(i),1)-s(trialbegind,1)>0));
      indtrialbegdel = [indtrialbegdel dummy];
      [m dummy] = max(find(s(rejectind(i),1)-s(winopenind,1)>0));
      indwinopendel = [indwinopendel dummy];
   end
   
   % the following are the: trialbeg and winopen column INDICES in s (not the times), corresponding to
   % good trials.
   indtrialbegkeep = setdiff(1:length(trialbegind),indtrialbegdel);
   indwinopenkeep = setdiff(1:length(winopenind),indwinopendel);
   trialbegindnew = trialbegind(indtrialbegkeep);
   winopenindnew = winopenind(indwinopenkeep);
else
   trialbegindnew=trialbegind;
   winopenindnew=winopenind;
end
stimbegindnew=zeros(size(trialbegindnew));
for h=1: max(size(trialbegindnew));
    stimbegindnew(h)=trialbegindnew(h)+find(s(trialbegindnew(h):wincloseind(h),2)==3232)-1;
end

%-------------------------------
% CUTTING DATAMATRIX INTO TRIALS
% BIG PROBLEM : acquisition (column LDind in s) and markers (column 2) indices
% of interest can be differents, because rex wrote the two datas separately !!
% We have to find data indices from marker indices by checking the time.
%-------------------------------
disp('removing zero data acquisition blocks')

% sdata contains the parts of s where laser doppler data [5th column] is not zero
sdataind = find(s(:,5));
% sdata timing should be monotonically increasing
f = find(s(sdataind(2:end),1)<s(sdataind(1:end-1),1));
if ~isempty(f)
    disp('(some data points go back in time)')
    badindices = [];
    f = [f ; length(sdataind)];
    for i=1:length(f)-1 % lets remove them
        pause(0)
        badindices = [badindices (f(i)+find(s(sdataind(f(i)+1:f(i+1)),1)<=s(sdataind(f(i)),1)))'];
    end
    pause(0)
    sdataind(badindices,:)=[];
end

% look for data indices
nstim = length(stimbegindnew);
stimbeginddata = ones(1,nstim);
nodatastim = [];
for i=1:nstim
    t = s(stimbegindnew(i),1);
    f = find_dicho(s(sdataind(:),1),t);
    if isempty(f) 
        disp(['attention: no data for stim ' num2str(i)])
        nodatastim(end+1) = i;
    else
        stimbeginddata(i) = f;
    end
end
% eliminating stims with no data
trialbegindnew(nodatastim) = [];
indtrialbegkeep(nodatastim) = [];
stimbegindnew(nodatastim) = [];
winopenindnew(nodatastim) = [];
wincloseindnew = wincloseind;
wincloseindnew(nodatastim) = [];
stimbeginddata(nodatastim) = [];

disp('cutting datamatrix into trials')
% calculating the length of the shortest trial (centered on stimbeg).
% Here, it is assumed that the time runs continously between daqwin
% opening and closing, thus delta(indices) are equivalent to delta(timesteps).
datapointsbefstimmin=min(s(stimbegindnew,1)-s(winopenindnew,1));
datapointsaftstimmin=min(s(wincloseindnew,1)-s(stimbegindnew,1));
% 1st dim of snew: data points, 2nd dim: the different columns of the B file, 3rd dim: the
%   different trials, without the rejected ones. Trial order: like stimbegindnew.
snew=zeros(datapointsbefstimmin+datapointsaftstimmin+1,size(s,2),length(stimbegindnew));
for i = 1:size(snew,3)
    snew(:,:,i)=s(sdataind(stimbeginddata(i)-datapointsbefstimmin:stimbeginddata(i)+datapointsaftstimmin),:);
end
%clear s

% SORTING TRIALS ACCORDING TO ECODES
disp('sorting trials according to Ecodes')
% find the indices of the different stimulusIDs in the trials: second dimension=stimID. first dimension= 
% trial index in trialbegindnew. The output is alrady sorted per stims as in
% stimlist (all non-rejected trials 1st conds, all trials 2nd cond. etc.)
stimindices=zeros(size(trialbegindnew,1),size(stimlist,2));
pause(0)
for j = 1:size(trialbegindnew,1)
    dummy=find(ismember(s(trialbegindnew(j):wincloseindnew(j),2),stimlist));
    if isempty(dummy), error('missing stimID-Ecodes'), end
    if length(dummy)>2, error('problem with stimID-Ecodes'), end
    i=find(stimlist==s(trialbegindnew(j)+dummy-1,2));
    stimindices(j,i)=s(trialbegindnew(j)+dummy-1,2);
end

% trnum contains the index (in trialbegindnew) of the trial of a given stimulus, and stimind contains the 
% index of that stimulus in stimlist!!!
[trnum,stind]=find(stimindices); % get rid of the zeros.
clear stimindices

%-------------
% LAST DETAILS
%-------------
disp('sorting, time axis, multiplicity')
% ordering the trials in the new datamatrix according to different
% stimuli. Ordering: as specified in stimlist.
conditions=snew(:,LDind,trnum); % keeping only the laser doppler (5th analog REX channel).

% setting time axis
timeaxis=squeeze(snew(:,1,trnum));
for i = 1:length(trnum)
    timeaxis(:,i)=timeaxis(:,i)-timeaxis(datapointsbefstimmin+1,i); % centering the timeaxis on stimbegin
end
clear snew;
timeaxis=mean(timeaxis,2);

% calculating the number of trials for each condition.
multiplicity=zeros(1,length(stimlist));
for i=1:size(stimlist,2)
    multiplicity(i) = length(find(stind==i));
end

%----------------------
% CHECK OBTAINED TRIALS
%----------------------
disp('check trials')
% colorations
% col = 'bgrcmk'; obj = 'ox+sdv^<>ph';
% DEC = [repmat(col',11,1) repmat(obj',6,1)];
DEC = ['cccggmmbbrryykkc' ; 'o+x+x+x+x+x+x+xs']';
% raw data
fn_imvalue
figure(6), clf
ispan = 1:ceil(length(s)/25000):length(s);
plot(s(ispan,1)/1000,s(ispan,3:end))
hold on
% markers
pause(0)
plot(s(trialbegind,1)/1000,zeros(size(trialbegind))+200,'*r')
plot(s(winopenind,1)/1000,zeros(size(winopenind)),'*c')
plot(s(stimbegind,1)/1000,zeros(size(stimbegind)),'*b')
plot(s(stimendind,1)/1000,zeros(size(stimendind)),'*b')
plot(s(wincloseind,1)/1000,zeros(size(wincloseind)),'*c')
plot(s(rejectind,1)/1000,zeros(size(rejectind)),'*k')
hp = zeros(1,length(stimlist));
for i=1:length(stimIDind)
    f = find(stimlist==s(stimIDind(i),2));
    hp(f) = plot(s(stimIDind(i),1)/1000,400,DEC(f,:),'linewidth',2);
end
try
    legend(hp,'blank1','100a','100b','200a','200b','300a','300b','600a','600b', ...
    '1000a','1000b','1200a','1200b','2000a','2000b','blank2')
end
% obtained trials do they fit raw data ?
i=0;
for k=1:length(multiplicity)
    for j=1:multiplicity(k)
        i=i+1;
        plot((s(sdataind(stimbeginddata(trnum(i))),1)+timeaxis(1:20:end))/1000, ...
            conditions(1:20:end,:,i),DEC(k,1))
        plot(s(sdataind(stimbeginddata(trnum(i))),1)/1000,conditions(1,:,i),'*b')
    end
end
drawnow

clear s sdata

%--------------------------------
function c = find_dicho(t,x)
% dichotomy search in a sorted vector

a=1; b=length(t)+1; c = floor((a+b)/2);
while c~=a
    if t(c)==x
        return
    elseif t(c)>x
        b=c;
    else
        a=c;
    end
    c = floor((a+b)/2);
end
if (c==1 || c==length(t)), c=[]; end
