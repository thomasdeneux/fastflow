function Y2 = fast_normalizeseqtest_min(Yname,ftav)
% function Y = fast_normalizeseqtest(Y,ftav)
% function Y = fast_normalize('Y',ftav)
%---
% achieves vessel registration (estimation of movement based on CS)
% and normalization (divide by previous frames. Parameter ftav specifies how many frames are averaged)
% in numerator and denominator

% passage par référence étrange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

fprintf('dividing:           ')
[nj ni nt] = size(Y);

Y2 = zeros(size(Y));

semi_ftav = floor(ftav/2);

for t=1:semi_ftav
    Y2(:,:,t)=min(Y(:,:,1:t+semi_ftav),[],3);
end;
for t=semi_ftav+1:nt-semi_ftav
    Y2(:,:,t)=min(Y(:,:,t-semi_ftav:t+semi_ftav),[],3);
end;
for t=nt-semi_ftav+1:nt
    Y2(:,:,t)=min(Y(:,:,t-semi_ftav:nt),[],3);
end;


% $$$ Y = reshape(Y,nj*ni,nt);
% $$$ for i=1:nt
% $$$     fprintf('\b\b\b\b\b\b\b\b\b\b%4i/%4i\n',i,nt)
% $$$         if i > nt-floor(ftav/2)
% $$$             Y(:,i)=mean(Y(:,[i-floor(ftav/2):i+floor(ftav/2)]-floor(ftav/2)),2);
% $$$         elseif i < floor(ftav/2)+1
% $$$             Y(:,i)=mean(Y(:,[i-floor(ftav/2):i+floor(ftav/2)]+floor(ftav/2)),2);
% $$$         else
% $$$             Y(:,i) = mean(Y(:,i-floor(ftav/2):i+floor(ftav/2)),2);
% $$$         end
% $$$ end
% $$$ Y = reshape(Y,nj,ni,nt);

