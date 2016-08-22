function s = fast_resdispselect(varargin)
% function s = fast_resdispselect([S[,s[,'correction']])
% function s = fast_resdispselect('default')

% default resdisp structure
s0 = struct( ...
    'data',     'result', ...
    'edge',     1, ...
    'algo',     'track', ...
    'trial',    1, ...
    'res',      'all', ...
    'oddeven',  'all', ...
    'points',   'all', ...
    'fzsub',    0, ...
    'smooth',   0);


% input
% (get default resdisp structure)
if nargin==1 && isequal(varargin{1},'default')
    s = s0; 
    return
end
% (fastflow object)
if nargin>=1
    S = varargin{1};
else
    S = struct;
    S.parameters.resultavgtrials = true;
    S.cond = [0 1 2];
end
% (initial structure)
if nargin>=2
    s = varargin{2};
else
    s = evalin('base','ans');
    if ~isstruct(s), s = s0; end
end
s = fn_structmerge(s0,s,'skip');

% structure correction
if nargin==3 && strcmp(varargin{3},'correction')
    s = correctstructure(S,s);
    return
end

% main function
s = userselect(S,s);

end

%--- correction ---
function s = correctstructure(S,s)

doavgtrials = S.parameters.resultavgtrials;
condlist = getcondlist(S);

% 'average' or 'timeavg'?
if doavgtrials
    badflag  = 'timeavg';
    goodflag = 'average';
else
    badflag  = 'average';
    goodflag = 'timeavg';
end
isalone = ~iscell(s.trial);
if isalone, s.trial = {s.trial}; end
for i=1:length(s.trial)
    if isequal(s.trial{i},badflag), s.trial{i} = goodflag; end
end
if isalone, s.trial = s.trial{1}; end

% cancel odd and even?
if ~doavgtrials
    s.oddeven = 'all';
end

% cancel stim and rest?
if ~ismember(s.res,condlist)
    s.res = condlist{1};
end

end

%--- core function ---
function s = userselect(S,s)

doavgtrials = S.parameters.resultavgtrials;
averageflag = fn_switch(doavgtrials,'average','timeavg');
condlist = getcondlist(S);
dostimrest  = length(condlist)>1;

% possible choices
choices = struct( ...
    'data',     {{'volume' 'volumef' 'result' 'lines' 'section' 'width' 'flow'}}, ...
    'edge',     {{'current' 'select'}}, ...
    'algo',     {{'current' 'track' 'radon' 'gabor'}}, ...
    'trial',    {{'current' averageflag 'select'}}, ...
    'res',      {condlist}, ...
    'oddeven',  {{'all' 'odd+even' 'odd' 'even'}}, ... different syntax to make it faster
    'points',   {{'image' 'all' 'middle' '2/5' '3/5' '4/5'}}, ...
    'fzsub',    {{0 5 10 20 50 100}}, ...
    'smooth',   {{0 1 3 10 30 100}});
if ~doavgtrials
    choices.oddeven = {'all'};
end

% graphic size parameters
h = 35; wspace = 2;
W = 500; H = 9*h;
currentrow = 1;
currentcol = 1+wspace;
wlabel   = 55; hlabel   = 15;
wcompare = 20; hcompare = 20;
wchoice  = 55; hchoice  = 25;
wnum     = 70; hnum     = 25;
wparset  =120; hparset  = 30;
wok      = 40; hok      = 25;
ncharmax = 20;

% open window
hf = figure(861);
clf(hf)
set(hf, ...
    'numbertitle','off','name','Select result display', ...
    'menubar','none','toolbar','none','resize','off') %, ...
    %'windowstyle','modal')
pos0 = get(0,'pointerlocation')-[W H]/2;
pos0 = max(pos0,[10 60]);
ss = get(0,'screenSize');
pos0 = min(pos0,ss(3:4)-[W+10 H+25]);
set(hf,'position',[pos0 W H])

% current state
% (do compare?)
F = fieldnames(s);
nfields = length(F);
docompare = s;
for i=1:nfields
    f = F{i};
    docompare.(f) = iscell(s.(f)); 
end
if ~docompare.edge, docompare.edge = (isnumeric(s.edge) && s.edge<0); end
docompare.oddeven = false;
% (which choices are currently selected?)
s2 = s;
num = struct('edge',[],'trial',[]);
num.trial = [];
doallparset = false;
for i=1:nfields
    f = F{i};
    x = s.(f); 
    if ~iscell(x), x = {x}; end
    if strcmp(f,'oddeven') && length(x)>1
        % special case for oddeven
        if ~isequal(x,{'odd' 'even'}), error merde, end
        x = {'odd+even'};
    end
    y = choices.(f);
    z = zeros(1,length(x));
    for j=1:length(x)
        ok = false;
        for k=1:length(y), if isequal(y{k},x{j}), ok=true; break, end, end
        if ok
            z(j) = k;
        elseif ismember(f,{'edge' 'trial'}) && isnumeric(x{j})
            z(j) = find(strcmp(y,'select'));
            num.edge  = union(num.(f),x{j});
        elseif strcmp(f,'algo') && any(findstr(x{j},'all'))
            doallparset  = true;
            docompare.(f) = true;
            x{j} = strrep(x{j},'all','');
            k = find(strcmp(x{j},y));
            z(j) = k;
        elseif strcmp(f,'trial')
            z(j) = find(strcmp(y,'current'));
        else
            error('value is not one of the possible choices')
        end
    end
    s2.(f) = z;
end                
% (display image?)
doimage = (isscalar(s2.points) ...
    && strcmp(choices.points(s2.points),'image'));
if doimage && any(struct2array(docompare))
    error('no multiple visualization in image mode')
end
% (time axis is time, no trial?)
dotimetrial = (~doavgtrials && isequal(s2.trial,'timeavg'));


% display menus
% (data)
buttonlabel data
hucompare.data = buttoncompare(1);
huchoices.data = buttonchoiceset(1);
newline
% (edge)
buttonlabel edge
hucompare.edge = buttoncompare(2);
huchoices.edge = buttonchoiceset(2);
hunum.edge = buttonnum(2);
newline
% (algo)
buttonlabel algo
hucompare.algo = buttoncompare(3);
huchoices.algo = buttonchoiceset(3);
putspace
huparset = buttonparset(3);
newline
% (trial)
buttonlabel trial
hucompare.trial = buttoncompare(4);
if dotimetrial, set(hucompare.trial,'enable','off'), end
huchoices.trial = buttonchoiceset(4);
hunum.trial = buttonnum(4);
newline
% (result type)
buttonlabel res
hucompare.res = buttoncompare(5);
huchoices.res = buttonchoiceset(5);
if ~dostimrest, set([hucompare.res huchoices.res],'enable','off'), end
newline
% (odd/even filtering)
buttonlabel oddeven
hucompare.oddeven = buttoncompare(6);
set(hucompare.oddeven,'enable','off')
huchoices.oddeven = buttonchoiceset(6);
if ~doavgtrials, set([hucompare.oddeven huchoices.oddeven],'enable','off'), end
newline
% (frame zero subtraction)
buttonlabel points
hucompare.points = buttoncompare(7);
huchoices.points = buttonchoiceset(7);
newline
% (frame zero subtraction)
buttonlabel fzsub
hucompare.fzsub = buttoncompare(8);
huchoices.fzsub = buttonchoiceset(8);
newline
% (frame zero subtraction)
buttonlabel smooth
hucompare.smooth = buttoncompare(9);
huchoices.smooth = buttonchoiceset(9);
% (buttonok)
ok = buttonok;

% terminate
waitfor(ok,'value',1)
if ~ishandle(hf), return, end
% (remove the special choices 'select')
k = find(strcmp(choices.edge,'select'));
s2.edge  = setdiff(s2.edge,k);
k = find(strcmp(choices.trial,'select'));
s2.trial = setdiff(s2.trial,k);
% (put the selected choices in s)
for i=1:nfields
    f = F{i};
    x = choices.(f)(s2.(f));
    % special cases
    switch f
        case 'algo'
            if doallparset
                % mark algos as 'alltrack', 'allgabor', 'allradon'
                for j=1:length(x)
                    if ~strcmp(x{j},'current'), x{j}=['all' x{j}]; end
                end
            end
        case {'edge' 'trial'}
            % add the selected edge or trial numbers
            x = [x num2cell(num.(f))]; %#ok<AGROW>
        case 'oddeven'
            if ~isscalar(x), error programming, end
            if strcmp(x,'odd+even'), x={'odd' 'even'}; end
    end
    if isscalar(x), x=x{1}; end
    % assign value in s
    s.(f) = x;
end
close(hf)

        
% BUTTONS DISPLAY
    function putspace()
        currentcol = currentcol+wspace;
    end

    function newline()
        currentrow = currentrow+1;
        currentcol = 1+wspace;
    end

    function u = button(wb,hb,varargin)
        u = uicontrol('parent',hf, ...
            'position',[currentcol H-currentrow*h+(h-hb)/2 wb hb], ...
            varargin{:});
        currentcol = currentcol+wb;
    end

    function buttonlabel(str)
        button(wlabel,hlabel,'style','text','string',str, ...
            'backgroundcolor',get(hf,'color'));
        putspace()
    end

    function u = buttoncompare(ii)
        ff = F{ii};
        u = button(wcompare,hcompare,'style','togglebutton', ...
            'string','M', ...
            'value',docompare.(ff),'enable',fn_switch(~doimage), ...
            'backgroundcolor',get(hf,'color'), ...
            'callback',@(hu,evt)setcompare(ii,get(hu,'value')));
        if nargin>=3
            set(u,'enable',fn_switch(doenable))
        end
        putspace()
    end

    function u = buttonchoiceset(ii)
        ff = F{ii};
        values = choices.(ff);
        nchoice = length(values);
        u = zeros(1,nchoice);
        for jj=1:nchoice
            if isnumeric(values{jj}), values{jj}=num2str(values{jj}); end
            if length(values{jj})>ncharmax, values{jj}=values{jj}(1:ncharmax); end
            u(jj) = button(wchoice,hchoice,'style','togglebutton', ...
                'string',values{jj}, ...
                'fontsize',10, ...
                'callback',@(hu,evt)setchoice(ii,jj,hu));
        end
        set(u(s2.(ff)),'value',1) 
    end

    function u = buttonnum(ii)
        ff = F{ii};
        numbers = num.(ff);
        u = button(wnum,hnum,'style','edit', ...
            'string',num2str(numbers,'%i '),'horizontalalignment','left', ...
            'backgroundcolor','w', ...
            'enable',fn_switch(~isempty(numbers)), ...
            'callback',@(hu,evt)setnum(ii,hu));
    end

    function u = buttonparset(ii)
        u = button(wparset,hparset,'style','listbox', ...
            'string',{'first parameters set only','all parameters sets'}, ...
            'fontsize',8, ...
            'value',1+doallparset,'enable',fn_switch(docompare.algo), ...
            'callback',@(hu,evnt)setparset(ii,logical(get(hu,'value')-1)));
    end

    function u = buttonok()
        currentcol = W - wok - 2*wspace;
        u = button(wok,hok,'style','togglebutton','string','OK');
    end

% CALLBACKS
    function setcompare(ii,b)
        ff = F{ii};
        docompare.(ff) = b;
        if ~b
            set(huchoices.(ff)(s2.(ff)(2:end)),'value',0)
            % disable selections when 'select' is unselected
            if any(strcmp(choices.(ff)(s2.(ff)(2:end)),'select'))
                set(hunum.(ff),'enable','off','string','')
            end
            s2.(ff) = s2.(ff)(1);
            % keep only first number when 'select' is selected
            if strcmp(choices.(ff)(s2.(ff)),'select') && isvector(num.(ff))
                num.(ff) = num.(ff)(1);
                set(hunum.(ff),'string',num2str(num.(ff)))
            end
        end
        % enable/disable selection of 'all parameters sets' for 'algo'
        if strcmp(ff,'algo')
            if b
                set(huparset,'enable','on')
            else
                set(huparset,'enable','off','value',1)
            end
        end
    end

    function setchoice(ii,jj,hu)
        ff = F{ii};
        b = get(hu,'value');
        % enable selection when 'select' is selected
        if strcmp(choices.(ff){jj},'select')
            if b
                num.(ff) = 1;
                set(hunum.(ff),'enable','on','string',1)
            else
                num.(ff) = [];
                set(hunum.(ff),'enable','off','string','')
            end
        end
        % 'image' is uncompatible with any multiple choice
        if strcmp(ff,'points')
            if strcmp(choices.(ff){jj},'image')
                doimage = true;
                % reset the 'points' line and disable multiple choices
                set(huchoices.points(s2.points),'value',0)
                set(hucompare.points,'value',0,'enable','off')
                docompare.points = false;
                % disable all multiple choices
                for kk=setdiff(1:nfields,ii)
                    g = F{kk};
                    set(hucompare.(g),'value',0,'enable','off')
                    setcompare(kk,false)
                end
            elseif doimage
                doimage = false;
                for kk=1:nfields
                    g = F{kk};
                    set(hucompare.(g),'enable','on')
                end
                set(hucompare.oddeven,'enable','off')
                if ~doavgtrials, set(hucompare.res,'enable','off'), end
                if dotimetrial, set(hucompare.trial,'enable','off'), end
            end
        end
        % 'all' might be uncompatible with other 'trial' choices
        if ~doavgtrials && strcmp(ff,'trial')
            if strcmp(choices.trial{jj},'timeavg')
                dotimetrial = true;
                set(hucompare.trial,'value',0,'enable','off')
                setcompare(ii,false)
            elseif dotimetrial
                dotimetrial = false;
                set(hucompare.trial,'enable',fn_switch(~doimage))
            end
        end
        % set/add/remove new choice
        if docompare.(ff)
            if b
                s2.(ff) = union(s2.(ff),jj);
            else
                s2.(ff) = setdiff(s2.(ff),jj);
            end
        else
            set(huchoices.(ff)(s2.(ff)),'value',0)
            % disable selections when 'select' is unselected
            if any(strcmp(choices.(ff)(s2.(ff)),'select'))
                num.(ff) = [];
                set(hunum.(ff),'enable','off','string','')
            end
            s2.(ff) = jj;
        end      
    end

    function setnum(ii,hu)
        ff = F{ii};
        str = get(hu,'string');
        numbers = str2num(str); %#ok<ST2NM>
        % check
        if (isempty(numbers) && ~isempty(str)) ...
                || (~docompare.(ff) && ~isscalar(numbers)) ...
                || any(mod(numbers,1)) ...
                || ((~strcmp(ff,'edge') || ~docompare.edge) && any(numbers<=0))
            % wrong value -> put the previous one
            set(hu,'string',num2str(num.(ff),'%i '))
            return
        end
        % save value
        num.(ff) = numbers;
    end

    function setparset(ii,b) %#ok<INUSL>
        doallparset = b;
    end

end

%------
% tools
%------

function condlist = getcondlist(S)

conds = unique(S.cond);
if isempty(conds) || isscalar(conds)
    condlist = {'all'};
    return
end
if ~ismember(0,conds), error('no blank condition'), end
stims = setdiff(conds,0);
condbase = {'all' 'rest' 'stim'};
if isscalar(stims)
    condmore = {};
    conddiv = {'stim/rest'};
else
    condmore = cell(1,length(stims));
    for i=1:length(stims)
        condmore{i} = ['stim' num2str(i)];
    end
    conddiv = {'stim/rest' sprintf('stim%i/stim%i',stims(1),stims(2))};
end
condlist = [condbase condmore conddiv];

end