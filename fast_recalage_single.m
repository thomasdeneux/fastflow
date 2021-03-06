function [ISO, Y] = fast_recalage_single(Yname,CS,iso,indices)
% same as fast_recalage, but heavy computations are made in class single,
% and Y output is of class single
%---
% function ISO = fast_recalage(Y,CS[,iso[,indices]])
% function [ISO Y] = fast_recalage('Y',CS[,iso[,indices]])
% function Y = fast_recalage('Y',ISO)
% function x2 = fast_recalage(x,ISO,nj,ni); <x(2,np),ISO(np,3)->x2(2,np,nt>
%---
% achieves vessel registration (estimation of movement based on CS)
% and normalization (divide by mean)

if nargin==0, help fast_recalage_single , return, end 

% image de points par la transformation
if ~ischar(Yname) && any(size(Yname)==2)
    x = Yname; if size(x,1)~=2, x=x'; end, x(3,:)=1;
    ISO = CS; if size(ISO,2)~=3, ISO=ISO'; end
    nj = iso;
    ni = indices;
    nt = size(ISO,1);
    np = size(x,2);
    x2 = zeros(2,np,nt);
    for i=1:nt
        M = fast_register_single(ISO(i,:),nj,ni);
        x2(:,:,i) = M*x;
    end
    ISO = x2;
    return
end

% passage par r�f�rence �trange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    if ~strcmp(class(Y),'single')
        error('data must be of class single when being passed as a reference to fas_recalage_single function')
    end
    evalin('caller',['clear ' Yname])
else
    Y = single(Yname);
end
[nj ni nt] = size(Y);

% resampling only
if nargin<3 && any(size(CS)==3)
    ISO=CS;
    if size(ISO,2)~=3, ISO=ISO'; end
    fn_progress('resampling',nt)
    for i=1:nt
        fn_progress(i)
        Y(:,:,i) = fast_register_single(Y(:,:,i),ISO(i,:));
    end
    ISO = Y;
    return
end

% coregister (and resample if second output required)
if nargin<3
    iso=[0 0 0];
end
if nargin<4
    skipe = 20;
    bin = 3;
    [I,J] = meshgrid(1+skipe:bin:ni-skipe,1+skipe:bin:nj-skipe);
    indices = J(:) + nj*(I(:)-1);
end
sigma = 2; % slight smoothing parameter

fn_progress('registering:',nt)
% register first frame to reference image using full window
fn_progress(1)
Y1 = filt2(Y(:,:,1),sigma);
iso = fast_register_single(filt2(CS,sigma),Y1,iso,indices);
[nj ni] = size(CS);
M1 = fast_register_single(iso,nj,ni);
% register all frames to the first frame using smaller window (defined in
% 'indices') ; this way, even if the registering to reference image is
% problematic, frames remain registered together 
ISO = zeros(nt,3); iso = [0 0 0];
exponent = 10^(floor(log10(nt))-2);
for i=1:nt
    if ~mod(i,exponent), fn_progress(i), end
    
    % register every frame to first frame
    iso = fast_register_single(Y1,filt2(Y(:,:,i),sigma),iso,indices);
    % deduce the registration of the frame to reference image
    M2 = fast_register_single(iso,nj,ni);
    ISO(i,:) = fast_register_single(M2*[M1;0 0 1],nj,ni);
    
    % resample if required
    if nargout>1, Y(:,:,i) = fast_register_single(Y(:,:,i),ISO(i,:)); end
end
fn_progress('end')


