function fast_visurecalage(ISO)
% function fast_visurecalage(ISO)
%---
% click on the traces to select bad trials (variable 'bad' in base
% workspace)

if nargin==0
    try
        ISO = evalin('base','ISO');
    catch
        help fast_visurecalage
    end
end

[nt dum nexp] = size(ISO);
iso = permute(ISO,[1 3 2]);
iso = reshape(iso,nt*nexp,3);

% display
figure(1)
fn_imvalue
taxis = (0:nt*nexp-1)/nt +.5;
hl=plot(taxis,[iso(:,1)*500 iso(:,2:3)]);
legend('rotation','translation x','translation y')
axis tight, ax = axis;
for i=0:nexp
    line([i i]+.5,ax(3:4),'color',[.5 .5 .5])
end

% selection of bad trials
assignin('base','bad',[]);
set(hl,'HitTest','on','ButtonDownFcn', [...
    'p=get(gca,''currentpoint'');' ...
    'k=round(p(1));' ...
    'bad = union(bad,k)'])



