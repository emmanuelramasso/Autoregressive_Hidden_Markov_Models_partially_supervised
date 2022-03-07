function v2 = duree_vs_amplitude_netb(mKmeans,X, T, Q, feat1, feat2)

global noms
global reordonnerPourVisu
if isempty(reordonnerPourVisu), reordonnerPourVisu = true; end

if size(mKmeans,1)==1 || size(mKmeans,2)==1
   bgkfcm = mKmeans;
else
   [u, bgkfcm]=max(mKmeans,[],2);
end


if reordonnerPourVisu %length(unique(v))==Q
   % REORDONNE
   v=bgkfcm;
   disp('Pour visu...')
   uu=1:Q;%unique(v);
   h=[];
   for j=1:length(uu)
      ff=find(v==uu(j));
      if length(ff)==0, ff=0; end
      h=[h ff(1)];% contient le premier instant d'apparition
   end
   clear ff
   [u w]=sort(h,'ascend');
   v2=zeros(size(v));
   for j=1:length(w)
      ff=find(v==w(j));
      v2(ff)=j;
   end
   bgkfcm=v2;
   clear v ff  u j h uu
else
    v2=bgkfcm;
end

C=zeros(length(T),Q);
%s={'r*','g*','k*','b*','c*','y*'};
% s={'r*','g*','b*','c*','k*','y*','m*','r*','b*','g*','r*','g*','b*','c*','k*','y*','m*','r*','b*','g*'};
% s={'*','s','o','x','d','v','^','<','>','p','h'};
% S=colormap('gray'); %close gcf
% f=(find(sum(S,2)>2.5)); S(f,:)=[];
% q=floor(size(S,1)/Q); q0=q; Z = zeros(Q,3);
% for i=1:Q, Z(i,:) = S(q0,:); Z(i+1,:)=S(floor(q0/2),:);q0=q0+q; end

l={};
% g2=figure; 
% set(gcf,'Position',[300 300 1200 700])
%scrsz = get(0,'ScreenSize')
figure; %set(gcf,'Position',[1 1 scrsz(3) 0.9*scrsz(4)])
Z=colormap('jet');
Z = Z(2:round(size(Z,1)/Q):end-1,:);
s={'*','s','o','x','d','v','^','<','>','p','h','*','s','o','x','d','v','^','<','>','p','h','*','s','o','x','d','v','^','<','>','p','h'};
 
% Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
NBL=2;
mgbas = 0.12; mggauche=0.06;
if mod(Q,2)==1, removeLast=true; else removeLast=false; end

ha = tight_subplot(round(Q/NBL),NBL,[.02 .05],[mgbas .03],[mggauche .03],removeLast)
   
for i=1:Q
   f=find(bgkfcm==i);
    %subplot(round(Q/2),2,i),
    axes(ha(i));  
    hold on, 
    set(gca,'FontName','Times','FontSize',28)
    
    prx = false; pry=false;
    if ismember(i,[1:NBL:Q]), pry=true; end
    if i==Q-1, prx = true; elseif i==Q && mod(Q,2)==0, pry=false; prx = true; end
    if ~((i==Q-1 && Q>2) || i==Q),  set(ha(i),'XTickLabel',''); end
    if ~pry, set(ha(i),'YTickLabel',''); end
    
    if nargin==4, 
       h=plot(X(:,13),log10(X(:,12)+1),'.','Markersize',10,'Linewidth',3)
       set(h,'Color',[0.85 0.85 0.85]);
       plot(X(f,13),log10(X(f,12)+1),s{i},'Color',Z(i,:),'Markersize',10,'Linewidth',3)%,'MarkerSize',12)
       %title('12 vs 13')
       if ~isempty(noms), if prx, xlabel(noms{13}),end, if pry, ylabel(noms{12}), end, 
       else  if prx, xlabel('Amplitude [dB]'), end, if pry, ylabel('Duration [µs]'); end, end
       set(gca,'Yscale','log')
    else
       %plot(X(:,feat1),X(:,feat2),'Color',[0.85 0.85 0.85],'LineStyle',':')
       h= plot(X(:,feat1),X(:,feat2),'.','Markersize',10,'Linewidth',3);
       set(h,'Color',[0.85 0.85 0.85]);
       plot(X(f,feat1),X(f,feat2),s{i},'Color',Z(i,:),'Markersize',10,'Linewidth',3)%,'MarkerSize',12)
       %title(sprintf('%d vs %d',feat1,feat2))
       if ~isempty(noms), if prx, xlabel(noms{feat1}), end, if pry, ylabel(noms{feat2}), end, 
       else  if prx, xlabel('Amplitude [dB]'), end, if pry, ylabel('Duration [µs]'); end, end
    end
    axis tight
    grid on

%     if nargin==4, 
%        xlabel('Amplitude (dB)')
%        ylabel('Duration (log)')
%     end
    % ENERGIE PAR CLUSTER
    %cs3=caracteristiques_logiciel(:,POS_ENERGIE);
    %x=1:length(cs3); x(f)=[]; cs3(x)=0;
    %figure(g3(i)), plot(T,cumsum(cs3),'m'), hold on, plot(T,cs2)
    %title(sprintf('cumul energie - GK fr Kmeans - C%d',i))
    l{end+1}=sprintf('C%d',i);
end
axis tight

