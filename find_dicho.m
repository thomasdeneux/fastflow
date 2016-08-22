%--------------------------------
function c=find_dicho(t,x)
% dichotomy search in a sorted vector

a=1; b=length(t)+1; c = floor((a+b)/2);
while c~=a
    if t(c)==x
        return
    elseif t(c)>x
        b=c;
    else
        a=c;
    end
    c = floor((a+b)/2);
end