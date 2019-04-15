%% bivariate KDE
x1 = Orte(:, 2);
x2 = Orte(:, 3);
% sampling distance in nm
sampDist = 100; % try different factor 100 at home sum(f)*sampDist * sampDist should be one
% grid spacing in nm
gridSpace = 50;
% bandwith for filter in nm
bandwith = 50;
xGridSpace = round((max(x1)-min(x1))/gridSpace);
yGridSpace = round((max(x2)-min(x2))/gridSpace);
[x1, x2] = meshgrid(min(Orte(:, 2)):sampDist:max(Orte(:, 2)), min(Orte(:, 3)):sampDist:max(Orte(:, 3)));
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
figure
[f, xi, bw] = ksdensity(Orte(:, 2:3), xi, 'Bandwidth', bandwith); % set density estimator to real data
% plotting
x = linspace(min(x1), max(x1), xGridSpace);
y = linspace(min(x2), max(x2), yGridSpace);
f = f * sampDist * sampDist; % f times sampling distance
[xq, yq] = meshgrid(x, y);
z = griddata(x1, x2, f, xq, yq); % transform z from nm into density per µm²
z = z ./ (max(x1)/(1000*gridSpace) * max(x2)/(1000*gridSpace)); % 1000 times gridSpace or y xGridSpace
imagesc(z.')
figure
h = contour(z.', 'ShowText', 'on');
set(gca, 'Ydir', 'reverse')
% plot histogram of f or z
[~,edges] = histcounts(log10(f));
figure
histogram(f,10.^edges)
set(gca, 'xscale','log')
% can be used for thresholding
% imagesc(z>0.00000000003)