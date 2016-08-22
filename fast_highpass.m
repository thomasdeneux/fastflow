function Y = fast_highpass(Yname,sigma)
% function Y = fast_highpass(Yname,sigma)


% passage par référence étrange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

h = fspecial('gaussian',ceil(1.5*sigma),sigma);

fn_progress('high pass:',size(Y,3))
for i=1:size(Y,3)
    fn_progress(i)
    %Y(:,:,i) = filt2(Y(:,:,i),sigma,'h');
    Y(:,:,i) = Y(:,:,i) - imfilter(Y(:,:,i),h);
    pause(0)
end

% [nj ni nt] = size(Y);
% if mod(sigma,1), error('sigma should be an integer'), end
% nj2 = ceil(nj/sigma);
% ni2 = ceil(ni/sigma);
% [subi subj] = meshgrid((1:ni2)*sigma-sigma/2,(1:nj2)*sigma-sigma/2);
% 
% disp('downsample')
% z = zeros(nj2,ni2,nt);
% for i=1:ni2
%     indi = min((1:sigma)+(i-1)*sigma,ni);
%     for j=1:nj2
%         indj = min((1:sigma)+(i-1)*sigma,nj);
%         z(j,i,:) = nanmean(reshape(Y(indj,indi,:),sigma*sigma,nt));
%     end
% end
% 
% disp('low-pass filter on downsampled data')
% h = fspecial('gaussian',3,1);
% z = imfilter(z,h,'symmetric');
% 
% fn_progress('making high-passed data',nt)
% for i=1:nt
%     fn_progress(i);
%     Y(:,:,i) = Y(:,:,i) - interp2(subi,subj,z,1:ni,1:nj);
% end
% 
