function addSpiceContours

ax = gca;
saLim = ax.XLim;
ctLim = ax.YLim;

sa = linspace(saLim(1),saLim(2));
ct = linspace(ctLim(1),ctLim(2));

[SA,CT] = meshgrid(sa,ct);
spice = gsw_spiciness0(SA,CT);

hold on;
[c,h] = contour( SA,CT,spice, 0:0.5:2, '--','color',0.4*[1,1,1]);
% clabel(c,h);

end