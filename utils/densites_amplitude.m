function Sortie = densites_amplitude(mkmeans,C4, leF, plotAREA, lexlabel)

if nargin==2, leF=13; plotAREA=1; lexlabel='Amplitude [dB]'; end
if nargin==3, plotAREA=1; lexlabel='Amplitude [dB]'; end
if nargin==4, lexlabel='Amplitude [dB]'; end
if size(mkmeans,2)==1 || size(mkmeans,1)==1
   bgkfcm=mkmeans; 
else [u, bgkfcm]=max(mkmeans,[],2);
end
u=unique(bgkfcm);
C=length(u);
% if plotAREA, 
%    s={'*','s','o','x','d','v','^','<','>','p','h'};
% %    S=colormap('gray'); close gcf
% %    f=(find(sum(S,2)>2.5)); S(f,:)=[];
% %    q=floor(size(S,1)/C); q0=q; Z = zeros(C,3);
% %    for i=1:C, Z(i,:) = S(q0,:); Z(i+1,:)=S(floor(q0/2),:);q0=q0+q; end
% %    B=Z; clear Z
% Z=colormap('jet');
% B = Z(2:round(size(Z,1)/Q):end,:);clear Z
 s={'*','s','o','x','d','v','^','<','>','p','h'};
% else
%    s={'rx','gs','bd','cp','k*','yv','m<','r>','b^','gh'};
%    
% end

l={};
scrsz = get(0,'ScreenSize')
figure; set(gcf,'Position',[1 1 scrsz(3) 0.9*scrsz(4)]), hold on
Z=colormap('jet');
Z = Z(2:round(size(Z,1)/C):end-1,:);

for i=1:C
   f=find(bgkfcm==u(i));
   if length(f)>0
      [x y]=ksdensity(C4(f,leF),...
         linspace(min(C4(:,leF)),max(C4(:,leF)),500),'width',3);% T(f));%'npoints',1000);%T(f));
      if nargout==1, Sortie{i} = [x ; y]; end
      if plotAREA, area(y,x,'FaceColor',Z(i,:)), alpha(0.8)
      else          plot(y,x,s{i}), end
      l{end+1}=sprintf('Cluster %d',u(i));
   end
end
h=legend(l,'Location','Best');
% set(h,'Position',[ 0.81954      0.34235      0.16814      0.37985])
% set(gcf,'Position',[1          26        1920         945])
set(gca,'FontName','Times','FontSize',28)
% set(gca,'Position',[0.081771         0.11      0.85729      0.85614])
xlabel(lexlabel,'FontName','Times','FontSize',28),
ylabel(sprintf('Probability'),'FontName','Times','FontSize',28)

axis tight
