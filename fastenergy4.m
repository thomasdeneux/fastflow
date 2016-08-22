function E = fastenergy4(dec,Y,t,x,twidth,xwidth)
% function E = fastenergy4(dec,Y,t,x[,twidth,xwidth])

% Input
if nargin<6
    twidth = 3;
    xwidth = 7;
end

% reject points too close to borders
[nt nx] = size(Y);
if t<=twidth || t>nt-twidth || x<=2*xwidth*twidth || x>=nx-2*xwidth*twidth
    E = 0;
    return
end

% make matrix A as a neighbourghood of (t,x) in Y, and shifting rows by dec

A = zeros(1+2*twidth,1+2*xwidth);
for i=-twidth:twidth
    q = floor(i*dec); r = i*dec-q;
    A(1+twidth+i,:) = (1-r)*Y(t+i,x+q+(-xwidth:xwidth)) ...
        + r*Y(t+i,x+q+1+(-xwidth:xwidth));
end
A = fn_add(A,-A(2,:));
E = sum(A(:).^2);


