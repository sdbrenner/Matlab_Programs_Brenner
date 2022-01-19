function tH = panelLabel(n,offsetX,offsetY)
% PANELLABEL creates a letter label (a), (b), (c), etc. in the upper left
% corner of a figure
%
%   panelLabel(n) creates a label with the n-th letter of the alphabet;
%   e.g., for n=1 the label will be (a), for n = 12 the label will be (l)
%
%   panelLabel(n,offset) allows for specification of the offset distance
%   from the top left corner of the figure (in inches).  If
%   offset is not specificied, a default value of 0.1 inches is used.
%
%   panelLabel(n,offsetX,offsetY) allows separate specification of both the
%   x- and y-offset.
%
%   tH = panelLable(...) returns the text handle associated with the label.

    if nargin < 2
        offsetX = 0.1;
    end
    if nargin < 3
        offsetY = offsetX;
    end

    labels = sprintfc('(%c)','a':'z');
    ax = gca;

    oldUnits = ax.Units;
    ax.Units = 'inches';
    posIn = ax.Position;
    tX = offsetX/posIn(3);
    tY = 1 - offsetY/posIn(4);

    tH = text( tX, tY, labels{n},....
          'units','normalized',...
          'HorizontalAlignment','left','VerticalAlignment','top',...3
          'FontName',ax.FontName,...
          'FontSize',ax.FontSize);
    ax.Units = oldUnits;
    
end