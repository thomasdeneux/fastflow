function fast_movie(Y)

CS = double(Y(:,:,1));
nt = size(Y,3);
hf = figure(794);
set(794,'numbertitle','off','name','MOVIE', ...
    'tag','fast movie'); % protect from fn_imvalue
cliptool = [2 (max(CS(:))+min(CS(:)))]/(max(CS(:))-min(CS(:)));
im = imagesc(CS'*cliptool(1) - cliptool(2),[-1 1]);
axis image
set(gca,'CLimMode','manual','xtick',[],'ytick',[])
hshift = uicontrol('style','togglebutton', ...
    'string','norm','value',1, ...
    'units','normalized','position',[.02 .015 .08 .05], ...
    'backgroundcolor',[.8 .8 .8]);
hclip = uicontrol('style','slider', ...
    'units','normalized','position',[.11 .015 .15 .05], ...
    'min',1,'max',2.5,'value',2);
hspeed = uicontrol('style','slider', ...
    'units','normalized','position',[.27 .015 .15 .05], ...
    'min',0,'max',3,'value',2);
hu = uicontrol('style','slider', ...
    'units','normalized','position',[.43 .015 .55 .05], ...
    'min',1,'max',nt,'value',1);
while ishandle(hf)
    speedfact = get(hspeed,'value');
    k = round(get(hu,'value'));
    if speedfact, k = 1+mod(k,nt); end
    set(hu,'value',k)
    % frame
    fr = double(Y(:,:,k));
    % normalize?
    donormalize = get(hshift,'value');
    if donormalize
        framavrl=5;
        idx = max(1,k-framavrl):min(nt,k+framavrl);
        blah = fr./mean(Y(:,:,idx),3) - 1;
        blah = blah * (10^get(hclip,'value'));
    else
        sigma = 4;
        blah = filt2(fr,sigma,'hm');
        blah = blah * (cliptool(1)*10^get(hclip,'value')/2);
    end
    % display
    set(im,'cdata',blah')
    pause(10^(-speedfact))
end
