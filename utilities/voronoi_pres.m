%{
bins = 0:0.00001:0.0002;
h1 = histcounts(cellDC.voronoi_micro_area, bins);
h2 = histcounts(cellPDC.voronoi_micro_area, bins);
h3 = histcounts(cellUDC.voronoi_micro_area, bins);
figure
h = bar(bins(1:end-1), [h1;h2;h3]);
%set(gca,'yscale','log')
legend('DC  #~180,000', 'PDC #~180,000', 'UDC #~180,000')
xlabel("Area in μm^{2}")
ylabel("Number of voronoi cells")
title({"Area histogram","Range 0 - 0.0002 μm^{2}"})
set(gca,'FontSize',18)
%}
f = figure;
[yvals,xvals] = ecdf(cellDC.voronoi_micro_densities);
semilogx(xvals,yvals,'-b')
hold on 
[yvals,xvals] = ecdf(cellPDC.voronoi_micro_densities);
semilogx(xvals,yvals,'-r')
hold on 
[yvals,xvals] = ecdf(cellUDC.voronoi_micro_densities);
semilogx(xvals,yvals,'-','color',[1.0, 0.7, 0])
title({"Cumulative density plot";})
xlim([100,1e6]);
ylabel("Cumulative fraction")
xlabel("Voronoi density in μm^{-2}")
lgd = legend('DC  #~180,000', 'PDC #~180,000', 'UDC #~180,000');
lgd.Location = 'northwest';
set(gca,'FontSize',18)


