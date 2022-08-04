function[]=printscript

% function to configure plotem figure for printing

c=get(gcf,'children');
c=sort(c);
for x=1:length(c)
	set(c(x),'Units','normalized')
	set(gcf,'PaperUnits','normalized','PaperPosition',[.04 .04 .92 .92])
end