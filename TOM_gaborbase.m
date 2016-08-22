function e=log_gabor_pyramid_wedge(s,periods,extents,ellips,angles)
% function e=log_gabor_pyramid_wedge(s,periods,extents,ellips,angles)
%---
% Compute a basis of Gabor filters (see
% http://en.wikipedia.org/wiki/Gabor_filter).
%
% g(x,y) = exp(-(u^2 + (v*GAMMA)^2)/SIGMA^2) * exp(2*pi*i*u/LAMBDA),
%
% where u = x*sin(THETA) - y*cos(THETA)
% and   v = x*cos(THETA) + y*sin(THETA).
%
% Input:
% - s           filter size
% - periods     vector of periodes (LAMBDA)
% - extents     vector of spatial attenuations (SIGMA)
% - ellips      vector of ellipticities (GAMMA, <1 to elongate in the
%               direction of the lines) 
% periods, extents and ellips should be all the same length, but it is also
% possible to give single values
% - angles      vector of directions (THETA), in radians!
%
% Output:
% - e           structure
%   e.Fk        s x s x (#periods*#scales) x #angles array: fft
%               representation of the filter basis

% develop
if nargin<1
    s = 512;
    periods = 4*power(1.5,[0 1 2 3]);
    extents = 4*power(2,[0 1 2 3]);
    angles = (30:5:70)*pi/180;
    ellips = [1 2];
end

% size
e.s=s;
e.R=e.s/2;

% Gabor parameters
n_period = length(periods);
n_extent = length(extents);
n_ellips = length(ellips);
e.n_scale = max([n_period n_extent n_ellips]);
if ~all(ismember([n_period n_extent n_ellips],[1 e.n_scale]))
    error('number of periods, extents and ellipcities should all be the same or be 1')
end
if n_period==1, periods = periods*ones(1,e.n_scale); end
if n_extent==1, extents = extents*ones(1,e.n_scale); end
if n_ellips==1, ellips = ellips*ones(1,e.n_scale); end
e.n_orient=length(angles);
e.M=(e.n_scale.*e.n_orient)*e.s*e.s; % overcompleteness
e.min_angle = angles(1);
e.max_angle = angles(end);

% coordinates
[e.fY,e.fX] = meshgrid(1:e.s,1:e.s);
[e.Y,e.X] = meshgrid(-(e.s)/2:(e.s)/2-1,-(e.s)/2:(e.s)/2-1);

% the filters
e.Gabor=zeros(e.s,e.s,e.n_scale,e.n_orient);
e.Fk=zeros(e.s,e.s,e.n_scale,e.n_orient);
fn_progress('building filter',e.n_scale*e.n_orient)
i_total = 0;
for i_orient=1:e.n_orient   
    theta = angles(i_orient);
    U = e.X*sin(theta) - e.Y*cos(theta);
    V = e.X*cos(theta) + e.Y*sin(theta);
    for i_scale=1:e.n_scale
        i_total = i_total+1;
        fn_progress(i_total)
        lambda = periods(i_scale);
        sigma = extents(i_scale);
        gamma = ellips(i_scale);
        ATTENUATION = exp(-(U.^2 + (V*gamma).^2) / sigma^2);
        PHASE = exp(2*pi*1i*U/lambda);
        GABOR = ATTENUATION .* PHASE;
        N_F=sum(sum(GABOR.*conj(GABOR)))/e.s/e.s;
        e.Gabor(:,:,i_scale,i_orient) = fftshift(GABOR)/sqrt(N_F);
        e.Fk(:,:,i_scale,i_orient) = fft2(fftshift(GABOR))/sqrt(N_F);
        % figure(101), imagesc(ATTENUATION)
        % figure(102), imagesc(imag(GABOR))
        % pause
    end
end

% close(101:102)

