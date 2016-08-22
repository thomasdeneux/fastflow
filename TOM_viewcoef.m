fn_4Dview(1,e.data','2d')
figure(1), colormap gray
fn_4Dview(2,permute(coef,[2 1 4 3]),'2dplot')

% for each point, each Gabor was the best?
c = coef;
[nt nx nsc ndir] = size(c);
n = nt*nx;
c = reshape(c,[nt*nx nsc*ndir]);
[dum idx] = max(c,[],2);
idx = reshape(idx,nt,nx);
isc = 1+mod((idx-1),nsc);
idir = 1+floor((idx-1)/nsc);
angmin = filter_struct.min_angle;
angmax = filter_struct.max_angle;
nbr_orient = filter_struct.n_orient;
ang = angmin + (0:(nbr_orient-1))/(nbr_orient-1) * (angmax-angmin);
angcos = cos(ang);
angsin = sin(ang);
[jj ii] = meshgrid(1:nx,1:nt);
ind = reshape(1:nt*nx,nt,nx);

% % display of directions (arrows)
% figure(1)
% hold on
% sub = 1;
% indsub = ind(1:sub:end,1:sub:end);
% dirx = angsin(idir)*1.2*sub;
% diry = angcos(idir)*1.2*sub;
% colsc = 'bgrcmyk';
% coldir = 'bgrcmykw';
% for k=1:ndir
%     f = find(idir(indsub)==k & idir(indsub)~=1);
%     col = coldir(1+mod(k-1,length(coldir)));
%     indsubk = indsub(f);
%     quiver(jj(indsubk),ii(indsubk),dirx(indsubk),diry(indsubk),0,col,'hittest','off')
% end
% hold off

% display of "interpreted" image
phase = zeros(nt,nx);
f = find(idir(ind)~=1);
phase(ind(f)) = Ic(ind(f) + n*(idx(ind(f))-1));
fn_4Dview(3,real(phase)','2d'), colormap gray
phase = zeros(nt,nx);

