function plotLines(values, labels, ttl, xlbl, ylbl, widths)
nlines = size(values, 1);
cmap = hsv(nlines);
lbls = cell(1,nlines);
figure; hold on;
for i=1:nlines
    plot(values(i,:), '-', 'Color', cmap(i,:), 'linewidth', widths(i));
    lbls{i} = [labels, ' ', num2str(i)];
end
legend(lbls);
title(ttl);xlabel(xlbl); ylabel(ylbl);
end