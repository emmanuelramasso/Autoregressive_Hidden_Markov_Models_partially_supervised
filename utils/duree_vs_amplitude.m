function v2 = duree_vs_amplitude(mKmeans,X, T, Q, feat1, feat2)

if nargin==4
   v2 = duree_vs_amplitude_netb(mKmeans,X, T, Q);
else
   v2 = duree_vs_amplitude_netb(mKmeans,X, T, Q, feat1, feat2);
end

return
% 
% if size(mKmeans,1)==1 || size(mKmeans,2)==1
%    bgkfcm = mKmeans;
% else
%    [u, bgkfcm]=max(mKmeans,[],2);
% end
% 
% 
% if 1%length(unique(v))==Q
%    % REORDONNE
%    v=bgkfcm;
%    disp('Pour visu...')
%    uu=1:Q;%unique(v);
%    h=[];
%    for j=1:length(uu)
%       ff=find(v==uu(j));
%       if length(ff)==0, ff=0; end
%       h=[h ff(1)];% contient le premier instant d'apparition
%    end
%    clear ff
%    [u w]=sort(h,'ascend');
%    v2=zeros(size(v));
%    for j=1:length(w)
%       ff=find(v==w(j));
%       v2(ff)=j;
%    end
%    bgkfcm=v2;
%    clear v ff  u j h uu
% end
% 
% C=zeros(length(T),Q);
% %s={'r*','g*','k*','b*','c*','y*'};
% s={'r*','g*','b*','c*','k*','y*','m*','r*','b*','g*','r*','g*','b*','c*','k*','y*','m*','r*','b*','g*'};
% l={};
% g2=figure; hold on,
% % set(gcf,'Position',[300 300 1200 700])
% scrsz = get(0,'ScreenSize')
% figure; set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])
% 
% for i=1:Q
%    f=find(bgkfcm==i);
%    subplot(round(Q/2),2,i), hold on,
%    set(gca,'FontName','Times','FontSize',24)
%    if nargin==4,
%       plot(X(:,13),log10(X(:,12)+1),'Color',[0.65 0.65 0.65],'LineStyle','.')
%       plot(X(f,13),log10(X(f,12)+1),s{i})%,'MarkerSize',12)
%       title('12 vs 13')
%       set(gca,'Yscale','log')
%    else
%       plot(X(:,feat1),X(:,feat2),'Color',[0.65 0.65 0.65],'LineStyle','.')
%       plot(X(f,feat1),X(f,feat2),s{i})%,'MarkerSize',12)
%       title(sprintf('%d vs %d',feat1,feat2))
%    end
%    axis tight
%    if nargin==4,
%       xlabel('Amplitude (dB)')
%       ylabel('Duration (log)')
%    end
%    % ENERGIE PAR CLUSTER
%    %cs3=caracteristiques_logiciel(:,POS_ENERGIE);
%    %x=1:length(cs3); x(f)=[]; cs3(x)=0;
%    %figure(g3(i)), plot(T,cumsum(cs3),'m'), hold on, plot(T,cs2)
%    %title(sprintf('cumul energie - GK fr Kmeans - C%d',i))
%    l{end+1}=sprintf('C%d',i);
% end
% axis tight

