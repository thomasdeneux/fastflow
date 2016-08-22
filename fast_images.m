fn_cd mars
load save/8nov_fastflow/divers
load save/8nov_fastflow/fastflow
load save/8nov_fastflow/doppler

% load save/8nov_fastflow/test
% Y = fast_loaddata(files(1,:),400,'normalize');

% ax1 = [176 208 232 264]; % t = 97 ou 196
% ax2 = [145  185  287  327]; % t = 196

% figure(1), clf, colormap gray, set(1,'doublebuffer','on')
% set(1,'position',[360   264   1000   250],'paperposition',[5 5 14 3.8])
% for j=1:4
%     ha(j)=axes('position',[.02+.25*(j-1) .08 .21 .80]);
% end
% for i=97
%     for j=1:4
%         axes(ha(j))
%         im = Y(:,:,i-5*(4-j));
%         imagesc(im-mean(im(:)),[-.01 .01])
%         axis(ax1)
%         set(gca,'XTick',[],'YTick',[])
%         title([num2str((5*(j-1))*5) 'ms'])
%     end
% end
% axes(ha(1)), title('{\bf rect. 1}   0ms')
% set(ha,'Clim',[-.005 .007])
% 
% figure(2), clf, colormap gray, set(1,'doublebuffer','on')
% set(2,'position',[360   264   1000   250],'paperposition',[5 5 14 3.8])
% for j=1:4
%     ha(j)=axes('position',[.02+.25*(j-1) .08 .21 .80]);
% end
% for i=196
%     for j=1:4
%         axes(ha(j))
%         im = Y(:,:,i-5*(4-j));
%         imagesc(im-mean(im(:)),[-.01 .01])
%         axis(ax1)
%         set(gca,'XTick',[],'YTick',[])
%         title([num2str((5*(j-1))*5) 'ms'])
%     end
% end
% axes(ha(1)), title('{\bf rect. 1}   0ms')
% set(ha,'Clim',[-.005 .007])
% 
% figure(3), clf, colormap gray, set(1,'doublebuffer','on')
% set(3,'position',[360   264   1000   250],'paperposition',[5 5 14 3.8])
% for j=1:4
%     ha(j)=axes('position',[.02+.25*(j-1) .08 .21 .80]);
% end
% for i=203
%     for j=1:4
%         axes(ha(j))
%         im = Y(:,:,i-5*(4-j));
%         imagesc(im-mean(im(:)),[-.01 .01])
%         axis(ax2)
%         set(gca,'XTick',[],'YTick',[])
%         title([num2str((5*(j-1))*5) 'ms'])
%     end
% end
% axes(ha(1)), title('{\bf rect. 2}   0ms')
% set(ha,'Clim',[-.005 .007])

% fn_cd tredac
% saveas(1,'img_movie1.fig','fig')
% saveas(1,'img_movie1.eps','psc2')
% saveas(2,'img_movie2.fig','fig')
% saveas(2,'img_movie2.eps','psc2')
% saveas(3,'img_movie3.fig','fig')
% saveas(3,'img_movie3.eps','psc2')

% figure(4), clf, colormap gray
% set(4,'position',[360   264   560   560],'paperposition',[5 5 8 8])
% imagesc(CS(92:350,5:250))
% line(ax1([1 2 2 1 1])-4,ax1([3 3 4 4 3])-91,'color','black')
% line(ax2([1 2 2 1 1])-4,ax2([3 3 4 4 3])-91,'color','black')
% line(ax1([1 2 2 1 1])-4,ax1([3 3 4 4 3])-91,'color','white','linestyle','--')
% line(ax2([1 2 2 1 1])-4,ax2([3 3 4 4 3])-91,'color','white','linestyle','--')
% text('string','1','position',ax1([1 3])-[4 6+91],'color','white')
% text('string','2','position',ax2([1 3])-[4 6+91],'color','white')
% set(gca,'XTick',[],'YTick',[])
% colors = {'k','k--';'m','m--';'c','c--';'r','r--';'g','g--';'b','b--';}; 
% ind = [12 10 18 21 22 29];
% for i=1:6
%     x = edges(25,ind(i)).points2(indix{ind(i)},:); 
%     line(x(:,1)-4,x(:,2)-91,'color',colors{i,1},'linewidth',2)
% end
% saveas(4,'img_structure1.eps','psc2')

% figure(5), clf, colormap gray
% set(5,'position',[360   264   560   560],'paperposition',[5 5 8 8])
% imagesc(Y(92:350,5:250,188),[.987 1.008])
% line(ax1([1 2 2 1 1])-4,ax1([3 3 4 4 3])-91,'color','black')
% line(ax2([1 2 2 1 1])-4,ax2([3 3 4 4 3])-91,'color','black')
% line(ax1([1 2 2 1 1])-4,ax1([3 3 4 4 3])-91,'color','white','linestyle','--')
% line(ax2([1 2 2 1 1])-4,ax2([3 3 4 4 3])-91,'color','white','linestyle','--')
% text('string','1','position',ax1([1 3])-[4 6+91],'color','white')
% text('string','2','position',ax2([1 3])-[4 6+91],'color','white')
% set(gca,'XTick',[],'YTick',[])
% saveas(5,'img_structure2.eps','psc2')


% figure(6), clf, colormap gray
% set(6,'position',[360   100   500   500],'paperposition',[5 5 20 18])
% for j=1:4
%     ha(j) = axes('position',[.12+.225*(j-1) .33 .185 .65]);
% end
% ha(5) = axes('position',[.12 .08 .86 .21]);
% for i=21 % which trial is best ?
%     e = fast_vt(18,i);
%     e = fast_solve2D(e,30,5,'A2f');
%     [nt0 nx0] = size(e.data);
%     nt = nt0; nx = nx0;
%     [I J] = meshgrid(1:nx0,1:nt0);
%     [XX YY] = meshgrid(21:nx,1:nt); sub = YY+nt0*(XX-1);
%     axes(ha(1))
%     imagesc(21:nx,taxis/1000,e.data(sub))
%     %xlabel('space(mm)')
%     set(gca,'XTick',[])
%     ylabel('time(s)')
%     axes(ha(2))
%     imagesc(filty2(e.data(sub),30,'hm'),[-.01 .01])
%     set(gca,'XTick',[],'YTick',[])
%     axes(ha(3))
%     imagesc(21:nx,1:nt,filty2(e.data(sub),30,'hm'),[-.01 .01])
%     [XX YY] = meshgrid(21:7:nx,1:20:nt); sub2 = YY+nt0*(XX-1);
%     hold on, quiver(I(sub2),J(sub2),-e.fu(sub2),e.fv(sub2),.6,'b'), hold off
%     set(gca,'XTick',[],'YTick',[])
%     axes(ha(4))
%     imagesc(e.flux(sub),[.25 .95])
%     set(gca,'XTick',[],'YTick',[])
%     axes(ha(5))
%     plot(taxis/1000,(fn_normalize(mean(e.flux(sub),2),'one')-1)*100)
%     axis tight; tax = axis; 
%     hold on
%     plot(taxis/1000,(fn_normalize(conditions(:,dopplernum(i)),'one')-1)*100,'r--')
%     plot(taxis/1000,(fn_normalize(mean(-e.data(sub),2),'one')-1)*1000,'g-.')
%     hold off
%     legend('estimated velocity','laser doppler flow','blood volume')
%     axis tight; ax = axis; axis([tax([1 2]) ax([3 4])])
%     ylabel('% signal change')
%     xlabel('time (s)')
%     % ylabel('estimated velocity')
% end
% fn_cd tredac
% saveas(6,'img_lines.eps','psc2')

% figure(7), clf, colormap gray
% set(7,'position',[360   100   600   280],'paperposition',[5 5 20 8.7])
% for k=1:6
% %     ha(k) = axes('position',[.5 .12+.18*(k-1) .47 .15]);
%     i = mod(k,2)+1; j=floor((k+1)/2);
%     ha(k) = axes('position',[.06+.47*(i-1) .12+.3*(j-1) .45 .28]);
% end
% % ha(6) = axes('position',[.01 .03 .48 .96]);
% % axes(ha(6))
% % imagesc(CS(92:350,5:250)), axis image, set(gca,'XTick',[],'YTick',[])
% ind = [12 10 18 21 22 29];
% indkeep = 25:770;
% colors = {'k','k--';'m','m--';'c','c--';'r','r--';'g','g--';'b','b--';}; 
% dopplerflux = fn_normalize(mean(conditions,2),'one');
% for i=1:6
% %     axes(ha(6))
% %     x = edges(25,ind(i)).points2(indix{ind(i)},:); line(x(:,1)-4,x(:,2)-91,'color',colors{i,1},'linewidth',2)
%     axes(ha(i))
%     flux = fn_normalize(flux0(:,ind(i))+2*flux1(:,ind(i)),'one');
%     plot(taxis(indkeep)/1000,(flux(indkeep)-1)*100,colors{i,1}, ...
%         taxis(indkeep)/1000,(dopplerflux(indkeep)-1)*100,'k--')
%     axis tight, ax=axis; axis([taxis(indkeep([1 end]))/1000 ax([3 4])])
%     if i>2, set(gca,'XTick',[]), 
%     else xlabel('time (s)'), 
%     end
%     if mod(i,2), set(gca,'YTick',[]), 
%     else if i==4, ylabel('% signal change'), end
%     end
% end
% % figure(7) % putain le recalage temporel !!!
% % flux = fn_normalize(mean(test1,2)+mean(test2,2),'one');
% % axes(ha(6)), plot(taxis(indkeep)/1000,flux(indkeep),'k', ...
% %     timeaxis/1000-.5,dopplerflux,'k--')
% % axis tight, ax=axis; axis([taxis(indkeep([1 end]))/1000 ax([3 4])])
% % set(gca,'YTick',[])
% fn_cd tredac
% saveas(7,'img_result.eps','psc2')

% REPONSE AU STIM: ON ABANDONNE

% figure(9), clf
% set(9,'position',[360   264   500   160],'paperposition',[5 5 25 8])
% mtest = mean(test(indkeep,[1 5 6 7 11 12 13 17 18 20 22 25 26 27 28 30 31]),2);
% subplot(1,2,1)
% plot(taxis(indkeep)/1000,mtest,'b',taxis(indkeep)/1000,filty(mflux,40,'lm'),'b')
% axis tight, ax = axis; axis([tax([1 2]) ax([3 4])]), hold on
% xlabel('time (s)')
% subplot(1,2,2)
% dopplertest = 
% plot(timeaxis/1000-.5,doppler,'k',timeaxis/1000-.5,filty(doppler,200,'lm'),'k')
% axis tight, ax = axis; axis([tax([1 2]) ax([3 4])])
% xlabel('time (s)')
% fn_cd tredac
% saveas(9,'img_estimatedflow.eps','psc2')

% FILT = 60;
% figure(10), clf
% set(10,'position',[360   514   300 210],'paperposition',[5 5 10 7])
% axes('position',[.1 .15 .8 .83])
% li=0;
% li=li+1; hli(li)=plot(taxis/1000,100*detrend(mtest,'constant'),'b');
% hold on
% plot(taxis/1000,filty(100*detrend(mtest,'constant'),FILT),'b','linewidth',1.5)
% li=li+1; hli(li)=plot(taxis/1000,100*detrend(dopplertest,'constant'),'r--');
% plot(taxis/1000,filty(100*detrend(dopplertest,'constant'),FILT),'r--','linewidth',1.5)
% li=li+1; hli(li)=plot(taxis/1000,-100*detrend(mvtest,'constant')*10,'g-.');
% plot(taxis/1000,-filty(100*detrend(mvtest,'constant'),FILT)*10,'g-.','linewidth',1.5)
% axis tight, hold off
% ylabel('% signal change'), xlabel('time (s)')
% legend(hli,'estimated velocity','laser doppler flow','blood volume')
% % axes('position',[.6 .15 . .83])
% % li=0;
% % li=li+1; hli(li)=plot(taxis/1000,100*detrend(mtestbis,'constant'),'b');
% % hold on
% % plot(taxis/1000,filty(100*detrend(mtestbis,'constant'),FILT),'b','linewidth',1.5)
% % li=li+1; hli(li)=plot(taxis/1000,100*detrend(dopplertestbis,'constant'),'r--');
% % plot(taxis/1000,filty(100*detrend(dopplertestbis,'constant'),FILT),'r--','linewidth',1.5)
% % li=li+1; hli(li)=plot(taxis/1000,-100*detrend(mvtestbis,'constant')*10,'g-.');
% % plot(taxis/1000,-filty(100*detrend(mvtestbis,'constant'),FILT)*10,'g-.','linewidth',1.5)
% % axis tight, hold off
% % set(gca,'YTick',[]), xlabel('time (s)')
% % legend(hli,'estimated velocity','laser doppler flow','blood volume')
% fn_cd tredac
% saveas(10,'img_test.eps','psc2')

