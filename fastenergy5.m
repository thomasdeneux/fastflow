function [F, DF] = fastenergy5(flux,ya,yb,alpha,beta)
% function [F DF] = fastenergy5(flux,ya,yb,alpha,beta)
%----
% Input:
% - flux        nx vector 
% - ya,yb       two nx vectors (volume at time t and t+1)

% correlations
nx = length(flux);
if length(ya)~=nx || length(yb)~=nx, error('dimension mismatch'), end
tt = (1:nx)'+flux(:);
yb2 = interp1(1:nx,yb,tt);
yb2(tt<1) = yb(1);
yb2(tt>nx) = yb(nx);
Fflow = yb2(:)-ya(:);

% smoothness
Fsmooth = diff(flux(:));

F = [alpha*Fflow ; beta*Fsmooth];

if nargout>1
    % correlations
    eps = 1e-3;
    ttbis = (1:nx)'+flux(:)+eps;
    yb2bis = interp1(1:nx,yb,ttbis);
    yb2bis(ttbis<1) = yb(1);
    yb2bis(ttbis>nx) = yb(nx);
    Dyb2 = (yb2bis-yb2)/eps;
    DFflow = spdiags(Dyb2,0,nx,nx);
    
    % smoothness
    DFsmooth = spdiags([ones(nx,1) -ones(nx,1)],[0 1],nx-1,nx);
    
    DF = [alpha*DFflow ; beta*DFsmooth];
end
