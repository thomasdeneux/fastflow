
% YOU CAN CHOOSE WHICH ACTION TO ASSOCIATE WITH THE 'ACTION BUTTON'

% %--- Replace algo and copy files ---
% f1 = M.e.flow(M.kflow);
% f = M.f;
% parest_manage(M,'replace')
% tmp = fn_input('trialsoff','1:450');
% trialsoff = eval(['[' tmp ']']);
% if isequal(trialsoff,1:450), return, end
% copyfiles(f,f1,trialsoff)

% %--- Show line images for all trials concatenated ---
% x = getdataf(M.f,1:450);
% [np nt nexp] = size(x);
% y = ones(np+1,nt,nexp,'single')*min(x(:)); y(1:np,:,:) = x;
% y = reshape(permute(y,[1 3 2]),(np+1)*nexp,nt);
% figure(1), imagesc(y)

%--- i3e_resultonevessel ---
i3e_resultonevessel2

% %--- Images axis image --- 
% axis(M.res(1).ha,'image')
% axis(M.res(4).ha,'image')
% axis(M.res(5).ha,'image')

% %--- Re-display vessels/fake parameters ---
% vessels_display(M)

% %--- Set current edge active/inactive ---
% i = M.kedge;
% S.edges(i).active = ~S.edges(i).active;
% if S.fake, vessels_display(M), else vsl_color(M,i), end


