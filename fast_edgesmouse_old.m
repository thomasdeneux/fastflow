function [edges, indices] = fast_edgesmouse(Y,varargin)
% function [edges indices] = fast_edgesmouse(Y[,CS][,clip][,edges]) -> user mode
% function [edges indices] = fast_edgesmouse(edges,nj,ni)           -> compute points2
%                                                                       
%---
% click in window to get edges
%
% Input:
% - Y       the volume movie
% - CS      static volume map
% - clip    clipping used for displaying the movie (default = [.99 1.01])
% - edges   structure to be appended
%
% Output:
% - edges   structure with fields
%           . points    points selected with mouse
%           . points2   resampling of the points (according to a
%                       curvilinear abscissae)
% - indices necessited indices in Y to interpolate data on edges
% Thomas, 29 oct 2005

if nargin==0, help fast_edgesmouse, return, end

disp('[attention, clip defini dans le fichier et non par un argument]')
disp('[attention, appel a fn_imvalue]')

fn_imvalue image

% Input
edges = struct([]);
CS = Y(:,:,1);
clip = [.99 1.01];
if ~isstruct(Y)
    mouseflag = true;
    [nj ni nt] = size(Y);
    for i=1:nargin-1
        a = varargin{i};
        if isstruct(a)
            edges = a;
        elseif length(a)==2
            clip = a;
        else
            CS = a;
        end    
    end
else
    mouseflag = false;
    edges = Y;
    [nj ni] = deal(varargin{:});
end

%if nargin<2
if mouseflag
    % help display
    disp('[for security, debugger will start on error]')
    dbstop if error
    disp('click in window to get edges')
    disp('press ''n'' for a new edge')
    disp('press ''a'' for enabling new points to current edge')
    disp('press ''p'' for a new point to current edge')
    disp('press ''m'' to mirror the direction of current edge')
    disp('press ''c'' to cancel changes to current edge')
    disp('press ''e'' to end current edge')
    disp('press ''d'' to delete current edge')
    disp('grap a point to move it')
    disp('grap a point and press ''i'' to insert a new one')
    disp('grap a point and press ''r'' to remove it')
    disp('grap a point and press ''b'' to break an edge into two components')
    disp('grap an end point, move it to the end point of another edge and press ''g'' to glue two edges together')
    disp('press ''*'' to set/unset markers')
    disp('press ''q'' to quit and save')

    % prepare display
    figure(1), clf, colormap gray, set(1,'doublebuffer','on'), %fn_imvalue end
    ha1 = gca;
    im1 = imagesc(Y(:,:,1)','parent',ha1,clip); axis image
    set(ha1,'climmode','manual')
    figure(2), clf, colormap gray, set(2,'doublebuffer','on')
    ha2 = gca;
    im2 = imagesc(CS','parent',ha2); axis image
    marker = {'none','*'};
    figure(3), clf, ha3=gca; colormap gray
    
    % draw edges
    nedges = length(edges);
    [lines1 lines2] = deal(zeros(1,nedges));
    for i=1:nedges
        x = edges(i).points(:,[2 1])';
        lines1(i)=line(x(1,:),x(2,:),'parent',ha1,'hittest','on','color','yellow','marker',marker{1},'linestyle',':');
        lines2(i)=line(x(1,:),x(2,:),'parent',ha2,'hittest','on','color','yellow','marker',marker{2});
        clear x
    end
    
    % set callbacks
    set([1 2],'WindowButtonDownFcn','setappdata(2,''buttonpressed'',true)')
    set([1 2],'WindowButtonUpFcn','setappdata(2,''buttonpressed'',false)')
    set([1 2],'windowbuttonmotionFcn','1;')
    set([1 2 3],'KeyPressFcn','1;','CurrentCharacter',' ') % trick to activate the CurrentCharacter buffer
    
    i=0;
    ok = true;               % perform looping
    addingpoint = false;     % i'm I adding new point to the end of current line?
    pointselectedflag=false; % i'm I dragging a point?
    curline=[];              % index of current line in the 'edges' structure, as well as in 'line1' and 'line2' handle collections
    curind=[];               % index of current point in current line
    curx=[];                 % points of currenct line
    tmplines=[];             % lines joining last point and mouse cursor
    while (ok)
        
        set(gcf,'currentcharacter',' ')
        
        pause(.01);
        delete(tmplines), tmplines=[];
        i = mod(i,nt)+1;
        set(im1,'CData',Y(:,:,i)')
        p = get(gca,'currentpoint'); if ~isempty(p), p=p([1 3]); end
        
        if getappdata(2,'buttonpressed')    % point dragging
                
            if ~pointselectedflag           % new line selection
                ho = gco;
                if isempty(ho) || ~strcmp(get(ho,'type'),'line'), setappdata(2,'buttonpressed',false), continue, end
                pointselectedflag = true;
                if ~isempty(curline)
                    set([lines1(curline) lines2(curline)],'color','yellow')
                    edges(curline).points = curx([2 1],:)';
                    curline = []; addingpoint = false;
                end
                curline=[find(lines1==ho) find(lines2==ho)];
                curx=get(ho,'xdata'); curx(2,:)=get(ho,'ydata');
                d=curx-repmat(p',1,size(curx,2)); d=d(1,:).^2+d(2,:).^2; 
                [dum curind]=min(d);
                set([lines1(curline) lines2(curline)],'color','red')
            end      
            
            switch get(gcf,'currentcharacter')
                case 'i'
                    curx = curx(:,[1:curind curind:end]);
                case 'r'
                    if size(curx,2)==1, continue, end
                    curx = curx(:,[1:curind-1 curind+1:end]);
                    if curind>size(curx,2), curind=size(curx,2); end
                    setappdata(2,'buttonpressed',false);
                case 'b'
                    if curind==1 || curind==size(curx,2), continue, end
                    nedges = nedges+1;
                    lines1(nedges) = line(curx(1,curind+1:end),curx(2,curind+1:end),'parent',ha1,'hittest','on','color','yellow','marker',marker{1},'linestyle',':');
                    lines2(nedges) = line(curx(1,curind+1:end),curx(2,curind+1:end),'parent',ha2,'hittest','on','color','yellow','marker',marker{2});
                    edges(nedges).points = curx([2 1],curind+1:end)'; 
                    curx = curx(:,1:curind-1);
                    setappdata(2,'buttonpressed',false);
                case 'g'
                    if curind~=1 && curind~=size(curx,2), continue, end
                    if curind==1, curx=curx(:,end:-1:1); end
                    test = zeros(2,2,length(lines1));
                    for k=setdiff(1:length(lines1),curline)
                        testa=get(lines1(k),'xdata'); testb=get(lines1(k),'ydata');
                        test(:,:,k) = [testa([1 end]) ; testb([1 end])];
                    end
                    d=test-repmat(p',[1 2 size(test,3)]); d=d(1,:).^2+d(2,:).^2;
                    [dum minind]=min(d);
                    glueline=floor((minind-1)/2)+1;
                    gluex=get(lines1(glueline),'xdata'); gluex(2,:)=get(lines1(glueline),'ydata');
                    switch mod(minind-1,2)
                        case 0  % start of segment
                            curx = [curx(:,1:end-1) gluex];
                        case 1  % end of segment
                            curx = [curx(:,1:end-1) gluex(:,end:-1:1)];
                    end
                    edges(glueline) = [];
                    nedges=nedges-1;
                    delete(lines1(glueline)), delete(lines2(glueline))
                    lines1(glueline)=[]; lines2(glueline)=[];
                    if glueline<curline, curline=curline-1; end
                    setappdata(2,'buttonpressed',false);                        
                otherwise % moving grabbed point
                    curx(:,curind) = p';
            end
            set(lines1(curline),'xdata',curx(1,:),'ydata',curx(2,:))
            set(lines2(curline),'xdata',curx(1,:),'ydata',curx(2,:))
            
        else
        
            if pointselectedflag
                pointselectedflag = false;
            else
                switch get(gcf,'currentcharacter')
                    case 'q'
                        if ~isempty(curline)
                            set([lines1(curline) lines2(curline)],'color','yellow')
                            edges(curline).points = curx([2 1],:)';
                            curline = []; addingpoint = false;
                        end
                        ok = false;
                    case 'a'
                        if isempty(curline), continue, end
                        addingpoint = true;
                    case 'n'
                        if ~isempty(curline)
                            set([lines1(curline) lines2(curline)],'color','yellow')
                            edges(curline).points = curx([2 1],:)';
                            curline = []; addingpoint = false;
                        end
                        nedges = nedges+1;
                        curline = nedges;
                        edges(curline).points = [];
                        curx = p';
                        lines1(curline) = line(curx(1,:),curx(2,:),'parent',ha1,'hittest','on','color','red','marker',marker{1},'linestyle',':');
                        lines2(curline) = line(curx(1,:),curx(2,:),'parent',ha2,'hittest','on','color','red','marker',marker{2});
                        addingpoint = true;
                    case 'p'
                        if isempty(curline), continue, end
                        curx(:,end+1) = p';
                        set(lines1(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                        set(lines2(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                    case 'd'
                        if isempty(curline), continue, end
                        delete(lines1(curline)), lines1(curline) = [];
                        delete(lines2(curline)), lines2(curline) = [];
                        edges(curline) = [];
                        nedges = nedges-1;
                        curline = []; addingpoint = false;
                    case 'm'
                        if isempty(curline), continue, end
                        curx = curx(:,end:-1:1);
                        set(lines1(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                        set(lines2(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                    case 'c'
                        if isempty(curline), continue, end
                        curx = edges(curline).points(:,[2 1])';
                        set(lines1(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                        set(lines2(curline),'xdata',curx(1,:),'ydata',curx(2,:))
                    case 'e'
                        if isempty(curline), continue, end
                        set([lines1(curline) lines2(curline)],'color','yellow')
                        edges(curline).points = curx([2 1],:)';
                        curline = []; addingpoint = false;
                    case '*'
                        if strcmp(marker{gcf},'*'), marker{gcf}='none'; else marker{gcf}='*'; end
                        set(lines1,'marker',marker{1})
                        set(lines2,'marker',marker{2})
                    otherwise
                        if ~addingpoint, continue, end
                        tmpx = [curx(:,end) p'];
                        tmplines(1)=line(tmpx(1,:),tmpx(2,:),'hittest','off','parent',ha1,'color','black');
                        tmplines(2)=line(tmpx(1,:),tmpx(2,:),'hittest','off','parent',ha1,'color','white','linestyle',':');
                        tmplines(3)=line(tmpx(1,:),tmpx(2,:),'hittest','off','parent',ha2,'color','black');
                        tmplines(4)=line(tmpx(1,:),tmpx(2,:),'hittest','off','parent',ha2,'color','white','linestyle',':');
                        continue
                end
            end
                        
            % display line figure
            if isempty(curline) || size(curx,2)<2, continue, end
            e = struct('points',curx([2 1],:)');
            [e.points2 L] = nicesnake(e.points,.8);
            e = fast_interpalongedges(e,Y);
            L = L*((length(e.points2)-1)/L(end))+1;
            imagesc(filtx(e.data,30,'hm')','parent',ha3)
            for i=1:length(e.points)
                line([.5 nt+.5],L([i i]),'color','k','parent',ha3)
            end

        end        
        
    end
    
    % nice snake
    for k=1:numel(edges)
        edges(k).points2 = nicesnake(edges(k).points,.8);
    end
    
end

    
% indices in Y
if nargout>1
    indices = [];
    for k=1:numel(edges)
        ind = floor(edges(k).points2);
        ind = ind(:,2) + nj*(ind(:,1)-1);
        ind = unique([ind ; ind+1 ; ind+nj ; ind+nj+1]);
        indices = [indices ; ind];
        if any(ind>ni*nj)
            keyboard
        end
    end
    indices = unique(indices);
end

%-------------------------
function [x2 L] = nicesnake(x,ds)

ni = size(x,1);
if ni==1, x2=x; return, end
L = zeros(ni-1,1);
for i=2:ni, L(i) = L(i-1)+norm(x(i,:)-x(i-1,:))+1e-3; end
if ~isempty(L), x2 = interp1(L,x,0:ds:L(end)); end
if norm(x(end,:)-x2(end,:))<ds/3
    x2(end,:)=x(end,:);
else
    x2 = [x2 ; x(end,:)];
end

%-------------------------
