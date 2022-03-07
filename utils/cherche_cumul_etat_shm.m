function [c w] = cherche_cumul_etat_shm(T,Q,v,tracerFig)


global PARAM_TRACE reordonnerPourVisu
if isempty(reordonnerPourVisu), reordonnerPourVisu = true; end
if isempty(PARAM_TRACE) || ~isfield(PARAM_TRACE,'xlabel'), lexlabel='Time (s)'; else lexlabel = PARAM_TRACE.xlabel;  end

%if nargin==4, error('NEW VERSION'), end

% s={'r','g','b','c','k','m','r','b','g'};
%s={'k:','k-','k--','k-.','k.'};
ok=reordonnerPourVisu;
if ~isempty(PARAM_TRACE)
   if isfield(PARAM_TRACE, 'reordonne')
      if PARAM_TRACE.reordonne==false
         ok=0;
      end
   end
end

if ok
   %if 1%length(unique(v))==Q
   % REORDONNE
   disp('Reordonne pour visu...')
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
   v=v2;
   clear v2 ff  u j h uu
else
   w = 1:Q;
end


% sequences de mode
%S = cherche_sequences_etats_et_confiance_2(v);

c=zeros(length(v),Q);
for i=1:length(v)
   if v(i)>0
      c(i,v(i))=1;
   end
end

if tracerFig>0
   if tracerFig==1, figure; end
   
   %scrsz = get(0,'ScreenSize');
   %scrsz = scrsz / 2
   %set(gcf,'Position',[1 1 scrsz(3) 0.9*scrsz(4)]);
   
   if isempty(PARAM_TRACE) || ~isfield(PARAM_TRACE,'cmap')
      col=colormap('jet');
      widthLine = 3;
      printCluster = true;
   else
      col=colormap(PARAM_TRACE.cmap);
      widthLine = PARAM_TRACE.width;
      printCluster = PARAM_TRACE.printCluster;
   end
   
   %    col=colormap('gray');
   col = col(2:round(size(col,1)/size(c,2)):end,:); save col col
   hold on,
   % k=1;
   l={};
   for i=1:size(c,2)
      cs= ((c(:,i)));
      %f=[find(diff(cs)~=0)]; % imprime ce qui change only
      %f=1:size(cs,1);
      cs = log10(cumsum(cs));     
      %g=find(~isinf(cs(f)));
      %hold on, plot([T(f(g(1))) ; T(f(g))],[0 ; cs(f(g))],'Color',col(i,:),'LineWidth',widthLine,'MarkerSize',5)
      hold on, plot(T,cs,'Color',col(i,:),'LineWidth',widthLine,'MarkerSize',5)
   end
   xlabel(lexlabel);%,'FontName','Times','FontSize',28),
   ylabel('log CSCA');%,'FontName','Times','FontSize',28)
   grid on
end
axis tight

%A=estimer_mattrans_fr_viterbi_0(Q, v);
%A = mk_stochastic(A)
%disp(s(1:Q))
