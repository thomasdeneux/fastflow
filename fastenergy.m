function F = fastenergy(FLUX,Yv0,Yt,CSU,CSV,alpha,beta)
% function F = fastenergy(FLUX,Yv0,Yt,alpha,beta)
%---
% Energie basée sur:
% - fit to data: flux optique (dI/dx v(x) + dI/dt = 0)
% - smoothness:  divergence du flux

Fflow = FLUX.*Yv0 + Yt;

Fsmooth = divergence(FLUX.*CSU,FLUX.*CSV);
% FXX = diff(FLUX.*CSU,1,2);
% FYY = diff(FLUX.*CSV,1,1);
% %Fsmooth = (FXX(1:end-1,:) + FXX(2:end,:) + FYY(:,1:end-1) + FYY(:,2:end) )/2;
% Fsmooth = [FXX(:) ; FYY(:)];

F = [alpha*Fflow(:) ; beta*Fsmooth(:)];


% % continuité / temps
% function F = fastenergy(FLUX,Yv0,Yt,alpha,beta)
% % function F = fastenergy(FLUX,Yv0,Yt,alpha,beta)
% 
% Fflow = FLUX.*Yv0 + Yt;
% Fsmooth = diff(FLUX);
% 
% F = [alpha*Fflow(:) ; beta*Fsmooth(:)];

