function flaglist1 = fast_flagcolors(flaglist)

if nargin<1, flaglist = evalin('base','S.flaglist'); end

nlist = length(flaglist.alllists);
list = flaglist.alllists(flaglist.k);
s = list.list;
nflag = length(s);

hline = 22;
curline = 1;
wtext = 25;     htext = 13;
wcolor = 60;    hcolor = 20;
wflag = 100;    hflag = 20;
wok = 25;       hok = 20;
wnewflag = 70;
wname = 90;
wchoice = 95;
wspace = 2;
curcolumn = 1;
W = wtext + wcolor + wflag + wspace*4;
H = (nflag+2)*(hline);

pos0 = get(0,'screenSize');
hf = figure(721);
figure(hf)
set(hf, ...
    'numbertitle','off','name','flags & colors', ...
    'menubar','none', ...
    'position',[pos0(3:4)/2-[W H]/2 W H]);
delete(findobj(hf,'type','uimenu'))

% Default output
flaglist1 = flaglist;

% Display buttons
namebutton = button('edit',wname,hok, ...
    'string',flaglist.alllists(flaglist.k).name,'backgroundcolor','w', ...
    'callback',@(u,evt)chgname(u));
choicebutton = button('popupmenu',wchoice,hok, ...
    'string',{flaglist.alllists.name '(new list)'},'value',flaglist.k, ...
    'callback',@(u,evt)chglist(u));
nextline
for i=1:nflag
    drawrow(i)
end
newflagbutton = button('pushbutton',wnewflag,hok,'string','new flag', ...
    'callback',@(u,evt)newflag());
okbutton = button('togglebutton',wok,hok,'string','ok');

% Wait for OK pressed
waitfor(okbutton,'value',1)
if ishandle(hf)
    close(hf)
    flaglist.alllists(flaglist.k).list = s;
    flaglist1 = flaglist;
end

% display sub-functions
    function redisplay()
        list = flaglist.alllists(flaglist.k);
        s = list.list;
        nflag = length(s);
        H = (nflag+2)*(hline);
        set(hf,'position',[pos0(3:4)/2-[W H]/2 W H]);
        allbuttons = findobj(hf,'type','uicontrol');
        rowbuttons = setdiff(allbuttons,[newflagbutton okbutton namebutton choicebutton]);
        delete(rowbuttons)
        for u = [namebutton choicebutton]
            pos = get(u,'pos');
            set(u,'pos',[pos(1) H-hline-(hok-hline)/2 pos(3) hok])
        end
        curline = 2; curcolumn = 1;
        for ii=1:nflag, drawrow(ii), end
        set(namebutton,'string',flaglist.alllists(flaglist.k).name)
        set(choicebutton,'string',{flaglist.alllists.name '(new list)'}, ...
            'value',flaglist.k)
    end

    function u = button(style,w,h,varargin)
        u = uicontrol('style',style, ...
            'position',[curcolumn H-curline*hline-(h-hline)/2 w h], ...
            varargin{:});
        curcolumn = curcolumn+w+wspace;
        if nargout==0, clear u, end
    end
    function nextline()
        curline = curline+1;
        curcolumn = 1;
    end
    function drawrow(iflag)
        button('text',wtext,htext, ...
            'string',num2str(iflag));
        button('edit',wcolor,hcolor, ...
            'string',col2str(s(iflag).color), ...
            'backgroundcolor',s(iflag).color, ...
            'callback',@(u,e)chgcolor(u,iflag))
        button('edit',wflag,hflag, ...
            'string',s(iflag).flag, ...
            'backgroundcolor','w','horizontalalignment','left', ...
            'callback',@(u,e)chgflag(u,iflag))
        nextline
    end

% action sub-functions
    function chgcolor(u,iflag)
        col = str2num(get(u,'string')); %#ok<ST2NM>
        if length(col)==3
            s(iflag).color = col;
            set(u,'backgroundcolor',col)
        end
        set(u,'string',col2str(s(iflag).color))
    end
    function chgflag(u,iflag)
        s(iflag).flag = get(u,'string');
    end
    function newflag()
        % new element to structure
        nflag = nflag+1;
        s(nflag).flag = '';
        s(nflag).color = floor(rand(1,3)*11)/10;
        flaglist.alllists(flaglist.k).list = s;
        % re-display
        redisplay()
    end
    function chgname(u)
        flaglist.alllists(flaglist.k).name = get(u,'string');
        redisplay();
    end
    function chglist(u)
        % save current list changes
        flaglist.alllists(flaglist.k).list = s;
        % switch to other list
        klist = get(u,'value');
        if klist==flaglist.k
            % nothing to do
        elseif klist<=nlist
            % switch to other (existing) list
            flaglist.k = klist;
            redisplay()
        else
            % create new list
            nlist = nlist+1;
            flaglist.k = nlist;
            flaglist.alllists(nlist) = struct('name','choose name', ...
                'list',struct('flag','','color', floor(rand(1,3)*11)/10,'defaultpars',[]));
            redisplay()
        end
    end
end



%---
function str = col2str(col)
    
    str = cell(1,3);
    for i=1:3
        if mod(col(i),.1)
            str{i} = ['.' num2str(round(100*col(i))) ' '];
        elseif mod(col(i),1)
            str{i} = ['.' num2str(10*col(i)) ' '];
        else
            str{i} = [num2str(col(i)) ' '];
        end
    end
    str = [str{:}];
    str(end)=[];
                
end