function [h ax1]=superpose_energie(E, T, boolCumsum, lexlabel)

if nargin==2, boolCumsum=true; lexlabel='Cumulated energy (aJ)'; end
if nargin==3, if boolCumsum, lexlabel='Cumulated energy (aJ)'; else lexlabel='Loading (kN)'; end, end

if 1
   ax1=gca;
   axis tight
   improve_figure
   
   if boolCumsum, s='m'; else s= 'k'; end
   
   ax2 = axes('Position', get(ax1,'Position'),...
      'YAxisLocation','right',...
      'Color','none',...
      'YColor',s,'FontName','Times','FontSize',28);
   
   if boolCumsum,   
      cs=cumsum(E);
   else
      cs=E;
   end
   
   %cs=cs/maxcs;
   %cs=cs*max(max(log10(cumsum(c))));
   
   if boolCumsum, s='m'; else s= [0.3 0.3 0.3]; end
   
   if boolCumsum
      hold on, h=plot(T, cs,'Color',s,'LineWidth',4,'Parent',ax2);
   else
      hold on, h=patch([T(:) ; nan], [cs(:) ; nan], s);
      set(h,'edgecolor',s,'LineWidth',4,'Parent',ax2,'edgealpha',0.4);
   end
   
   %axis tight
   ylabel(lexlabel,'FontName','Times','FontSize',28);
  
   % pas d'abcisse pour axe2
   %set(ax2,'XTick',get(ax1,'XTick'));
   set(ax2,'XTick',[]);
   %set(ax2,'Position', get(ax1,'Position'));
   %set(ax1,'Position', get(ax2,'Position'));
  
   % important pour avoir meme dilatation en x et y sur les 2 axes
   xlimits = get(ax1,'XLim');
   set(ax2,'XLim',xlimits);
  
  %get(ax2,'YLim')
   
   
   %g=get(ax1);
   %axis tight
   
else
   cs=cumsum(E);
   cs=cs/max(cs);
   cs=cs*max(max(log10(cumsum(c))));
   hold on, h=plot(T, cs,'m','LineWidth',4);
end


