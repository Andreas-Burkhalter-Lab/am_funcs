function hletters = LetterFigures(axpos,hfig,offset)
if (nargin < 2)
        hfig = gcf;
end
if (nargin < 3)
  offset = [-0.02 0.05];
end
hcfig = gcf;
hcax = gca;
figure(hfig);
hax = axes('Position',[0 0 1 1],'Visible','off');
nletters = length(axpos);
for i = 1:nletters
        tl(i,:) = [axpos{i}(1),axpos{i}(2)+axpos{i}(4)];
        curlet = char(double('A') + i-1);
        cpos = tl(i,:) + offset;
        hletters(i) = text(cpos(1),cpos(2),curlet,...
                'Visible','on','FontWeight','bold',...
                'FontSize',14);
end
% capital letters, bold, 14pt font
