function raisexlabel(pix,hfig)
% raisexlabel: shifts all xlabels vertically by given # of pixels
% raisexlabel(pix) or raisexlabel(pix,figurehandle)
if (nargin == 1)
        hfig = gcf;
end
hax = findobj(hfig,'Type','axes');
hxlabel = get(hax,'XLabel');
nlabel = length(hxlabel);
for i = 1:nlabel
        hl = hxlabel{i};
        utext = get(hl,'Units');
        set(hl,'Units','pixels');
        pos = get(hl,'Position');
        pos(2) = pos(2)+pix;
        set(hl,'Position',pos);
        set(hl,'Units',utext);
end
