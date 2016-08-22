

if ~exist('M','var')
    msgbox('Variable ''M'' not in base workspace, please click MADFLOW > ''object in basespace'' to enable access to results', ...
        'Warning','modal')
    return
end

hf = gcf;

%==========================================================================
% START USER EDIT HERE
%
% YOU CAN CHOOSE WHICH ACTION TO ASSOCIATE WITH THE 'USER DISPLAY'
%
% Useful tips:
% M.e is the current vessel
% M.f is the current 
%==========================================================================

x = M.u.tmpdata.section_current.data;
[np nt] = size(x);
outidx = [1:floor(np/4) ceil(3*np/4):np];
x = fn_mult(x,1./mean(x(outidx,:),1));

profile = fast_vesselwidth(x);
width = profile(1,:);
center = profile(2,:);

figure(1)
imagesc(x)
line([1:nt NaN 1:nt],[center+width/2 NaN center-width/2],'color','b')

figure(2)
plot(width)

%==========================================================================
% END USER EDIT
%==========================================================================

figure(hf)
