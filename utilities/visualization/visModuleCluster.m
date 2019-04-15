function [varargout] = visModuleCluster(locTable, visMethod, pixelSize)
%visModuleCluster visualization for cluster algorithm
%
% Following visualization methods are supported
% histogramBinning      [2D only]
% gaussianBlur          [2D only]
%
% input parameters
%   locTable        Matrix that contains point coordinates. Rows should
%                   correspond to detected signals. Columns should be
%                   arranged in the follwing order (column 1 - x coordiante, column 2 - y coordiante, , column 4 - x loc. prec., column 5 - y loc. prec.)
%   visMethod       any of the specified visualization methods as string
%                   for example 'histogramBinning'
%   pixelSize       final super-resolution pixel size
%
% requires Image Processing Toolbox, Statistics and Machine Learning
% Toolbox
%
%   by Jan Neumann, IMB Mainz, 12.02.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 100; % needed in gaussian blur if individual blurring method is used, specifies binning for loc. prec.

if nargin < 3
    multiWaitbar('CloseAll');
    error('Please specify at least visualization method and pixel size!');
end
%% start with visualization method
switch visMethod
    case 'histogramBinning'
        x1 = double(ceil(locTable(:, 1)./pixelSize));
        y1 = double(ceil(locTable(:, 2)./pixelSize));
        % generate histogram binned image
        SRimage = sparse(x1-min(x1)+1, y1-min(y1)+1, 1);
        SRimage = full(SRimage);
    case 'gaussianBlur'
        minPrecision = floor(min(min(locTable(:, 4:5))));
        maxPrecision = ceil(max(max(locTable(:, 4:5))));
        % get dimension of image
        x1 = double(round(locTable(:, 1)./pixelSize));
        y1 = double(round(locTable(:, 2)./pixelSize));
        % generate histogram binned image
        [xSize, ySize] = size(full(sparse(x1-min(x1)+1, y1-min(y1)+1, 1)));
        minX = min(x1);
        minY = min(y1);
        SRimage = zeros(xSize, ySize, 'single');
        prevStep = minPrecision;
        steps = linspace(minPrecision, maxPrecision, nSteps);
        steps(1) = [];
        for ii = steps
            % find points with specific loc prec.
            idx = (mean(locTable(:, 4:5), 2) <= ii) & (mean(locTable(:, 4:5), 2) > prevStep);
            tempPoints = locTable(idx, :);
            globalBlur = ((0.5*(ii - prevStep)) + prevStep) / pixelSize;
            % generate histogram binned with current selection of points
            x1 = double(round(tempPoints(:, 1)./pixelSize));
            y1 = double(round(tempPoints(:, 2)./pixelSize));
            tempImage = sparse(x1-minX+1, y1-minY+1, 1, xSize, ySize); % check if x and y are correct
            tempImage = full(tempImage);
            % images should have fixed size to be abble to be added
            % check version of matlab
            if verLessThan('matlab', '8.5')
                h = fspecial('gaussian', 2*ceil(2*globalBlur) + 1, globalBlur);
                SRimage = SRimage + imfilter(tempImage, h, 'conv');
                % alternative implement own filter
                % maskSize = 2*ceil(2*globalBlur) + 1;
                % ind = -floor(maskSize/2) : floor(maskSize/2);
                % [X, Y] = meshgrid(ind, ind);
                % h = exp(-(X.^2 + Y.^2)/(2*globalBlur*globalBlur));
                % h = h / sum(h(:));
                % SRimage = SRimage + conv2(tempImage ,h);
            else
                SRimage = SRimage + imgaussfilt(tempImage, globalBlur);
            end
            prevStep = ii;
        end
end
varargout{1} = SRimage;
end