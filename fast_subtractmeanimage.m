function [Y immean] = fast_subtractmeanimage(Y)

nt = size(Y,3);
immean = zeros(nt,1);

fn_progress('substracting mean',nt)
for t=1:nt
    fn_progress(t)
    y = Y(20:end-20,20:end-20,t);
    immean(t) = nanmean(y(:));
    Y(:,:,t) = Y(:,:,t)-immean(t);
end
