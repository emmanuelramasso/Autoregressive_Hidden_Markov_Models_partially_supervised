function hdle = tracer_energie_dans_cluster(clusters, energie, T)

cc=tabulate(clusters);
u=cc(find(cc(:,2)>0),1);

f=figure;
Z=colormap('jet'); close(f)
Z = Z(2:round(size(Z,1)/max(u)):end,:);
s={'*-','s-','o-','x-','d-','v-','^-','<-','>-','p-','h-'};
% scrsz = get(0,'ScreenSize');
hdle=figure; hold on, %set(hdle,'visible','off')
% set(gcf,'Position',[1 1 scrsz(3) 0.9*scrsz(4)])
grid on

r=max(1,ceil(sqrt(length(u))));
% l=ceil(length(u)/r);
c=(cumsum(energie)); 
hold on, plot(T,c,'m','LineWidth',4)
r={'Tot. cumulated energy'};

E=zeros(length(clusters),1);
for i=1:length(u)
   f = find(clusters==u(i));
   E(:)=0;
   E(f) = energie(f);
   
   plot(T, (cumsum(E)),s{u(i)},'Color',Z(u(i),:),'LineWidth',2,'MarkerSize',5)
      
   r{end+1}=sprintf('Cum. energy in clust. %d',i);
end

legend(r,'Location','Best');

%xlim([0 T(end)]);

% set(gca,'FontName','Times','FontSize',28)
% set(gca,'Position',[0.081771         0.11      0.85729      0.85614])
% xlabel('Time (s)','FontName','Times','FontSize',28),
xlabel('Time (s)');%,'FontName','Times','FontSize',28),
% ylabel(sprintf('Cumulated Energy'),'FontName','Times','FontSize',28)
ylabel(sprintf('Cumulated Energy'));%,'FontName','Times','FontSize',28)
improve_figure


return

% s={'r','g','b','c','k','y','m','r','b','g'};
% u=unique(clusters);
% figure, set(gcf,'Position',[8         213        1891         685]);
% r=max(1,ceil(sqrt(length(u))));
% l=ceil(length(u)/r);
% c=cumsum(energie);
% for i=1:length(u)
%    f=find(clusters==i);
%    subplot(r,l,i),plot(T(f), cumsum(energie(f)),s{i},'LineWidth',2)
%    hold on, plot(T,c,'m','LineWidth',2)
%    legend('Total cum. Energy',sprintf('Cum. Energy in Cluster %d',i))
% end
%
% set(gcf,'Position',[1          26        1920         945])
% set(gca,'FontName','Times','FontSize',28)
% % set(gca,'Position',[0.081771         0.11      0.85729      0.85614])
% xlabel('Amplitude range','FontName','Times','FontSize',28),
% ylabel(sprintf('Probability'),'FontName','Times','FontSize',28)

global USELOG
if isempty(USELOG), USELOG=0; end
%USELOG=1;

%s={'rx','gs','bd','cp','k*','m<','r>','b^','gh'};

%u=unique(clusters);
cc=tabulate(clusters);
u=cc(find(cc(:,2)>0),1);
% figure, set(gcf,'Position',[8         213        1891         685]);

Z=colormap('jet'); close(gcf)
Z = Z(2:round(size(Z,1)/max(u)):end,:);
s={'*-','s-','o-','x-','d-','v-','^-','<-','>-','p-','h-'};
scrsz = get(0,'ScreenSize');
hdle=figure; hold on, set(hdle,'visible','off')
set(gcf,'Position',[1 1 scrsz(3) 0.9*scrsz(4)])
grid on

r=max(1,ceil(sqrt(length(u))));
% l=ceil(length(u)/r);
if USELOG, c=log10(cumsum(energie));
else c=(cumsum(energie)); end
hold on, plot(T,c,'m','LineWidth',4)
r={'Tot. cumulated energy'};
%r={'Energie cumulée totale'};
if nargout==2
   energieDansChaqueCluster = zeros(length(energie),length(u));
   %energieDansChaqueCluster(1,:)=0;
end

E=zeros(length(energie),1);
for i=1:length(u)
   f=find(clusters==u(i));
   E(:)=0;
   E(f) = energie(f);
   
   if USELOG, plot(T, log10(cumsum(E)),s{u(i)},'Color',Z(u(i),:),'LineWidth',2,'MarkerSize',5)
   else plot(T, (cumsum(E)),s{u(i)},'Color',Z(u(i),:),'LineWidth',2,'MarkerSize',5)
   end
      
%    if USELOG, plot(T(f), log10(cumsum(energie(f))),s{u(i)},'Color',Z(u(i),:),'LineWidth',2,'MarkerSize',5)
%    else plot(T(f), (cumsum(energie(f))),s{u(i)},'Color',Z(u(i),:),'LineWidth',2,'MarkerSize',5)
%    end
   if nargout==2
      %energieDansChaqueCluster(f,u(i)) = cumsum(energie(f));
      energieDansChaqueCluster(:,u(i)) = cumsum(E);
   end
   r{end+1}=sprintf('Cum. energy in clust. %d',i);
   %r{end+1}=sprintf('Energie cum. ds clust. %d',u(i));
end
% if nargout==1
%    for i=2:size(energieDansChaqueCluster,1)
%       for j=1:size(energieDansChaqueCluster,2)
%          if isnan(energieDansChaqueCluster(i,u(j)))
%             energieDansChaqueCluster(i,u(j)) = energieDansChaqueCluster(i-1,u(j));
%          end
%       end
%    end
% end
% if USELOG & nargout==1
%       energieDansChaqueCluster = log10(energieDansChaqueCluster);
% end
h=legend(r,'Location','Best');
xlim([0 T(end)]);
% set(gcf,'Position',[1          26        1920         945])
% scrsz = get(0,'ScreenSize')
% figure; set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])

set(gca,'FontName','Times','FontSize',28)
% set(gca,'Position',[0.081771         0.11      0.85729      0.85614])
% xlabel('Time (s)','FontName','Times','FontSize',28),
xlabel('Time (s)','FontName','Times','FontSize',28),
% ylabel(sprintf('Cumulated Energy'),'FontName','Times','FontSize',28)
ylabel(sprintf('Cumulated Energy'),'FontName','Times','FontSize',28)
