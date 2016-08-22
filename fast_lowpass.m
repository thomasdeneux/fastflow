function Y = fast_lowpass(Yname,sigma)
% function Y = fast_lowpass(Yname,sigma)
%---
% spatial smoothing using a gaussian


% passage par r�f�rence �trange -> save memory !!
if ischar(Yname)
    Y = evalin('caller',Yname);
    evalin('caller',['clear ' Yname])
else
    Y = Yname;
end

h = fspecial('gaussian',min(ceil(3*sigma),5),sigma);

fn_progress('low pass:',size(Y,3))
for i=1:size(Y,3)
    fn_progress(i)
    Y(:,:,i) = imfilter(Y(:,:,i),h);
    pause(0)
end

