function [] = compactionParameter(obj, varargin)
%compactionParameter determines chromatin compaction of nuclei in
%                    localization microscopy images
%
% The algorithm is based on Irianto et al.: "Quantification of chromatin
% condensation level by image processing", Med. Eng. Phys., (2014).
%
%   input:
%       varargin    struct containing the following fields in the same order:
%                   isRandom(default = 0): when true use random data
%                   nSize(default = 15): bin size for image (in nm)
%                   nFilters (default = 6): apply n times the filterType
%                   filterType default('average'): type of smothing filter for nucleus extraction
%                   PixRedFactor (default = 4): image reduction factor
%                   SobelThresh (default = 0.09)
%                   maxExternalObject (default = 50000): total number of pixels of connected object that should be removed (in px)
%                   maxExternalObjectReduced (default = 3125): maximum size of objects that should be removed after perimeter removal, /16 due to resizing of image
%                   nPerimeter (default = 2): perimeter size of nucleus that will be removed
%                   plotOverlay (default= false)
%   output:
%                   (1) as m x 4 matrix, where m corresponds to a nuclei
%                   column 1: area
%                   column 2: circularity between 0 (circle) and 1 (line segment)
%                   column 3: number of edges
%                   column 4: chromatin compaction parameter
%
% by Amine Gourram (see master thesis)
% modified by Jan Neumann, MPI for Chemistry, 08.11.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init parameters
if nargin > 2
    error('Too many input arguments!')
end
if obj.dimension == 3
	error('Compaction parameter analysis is implemented for 2D data only');
end
optargs = {0, 15, 6, 'average', 4, 0.09, 50000, 3125, 2, false};
if ~isempty(varargin)
    numvarargs = numel(fieldnames(varargin{1}));
    fields = fieldnames(varargin{1});
    for ii = 1:numvarargs
        optargs{ii} = varargin{1}.(fields{ii});
    end
end
[isRandom, binSize, nFilters, filterType, PixRedFactor, SobelThresh, maxExternalObject,...
    maxExternalObjectReduced, nPerimeter, plotOverlay] = optargs{:};
if isRandom == true
    positions = obj.randomTable;
    dataMat = obj.randomClusterStruct;
else
    positions = obj.positionTable;
    dataMat = obj.clusterStruct;
end
multiWaitbar('Calculating CCP. Please wait!', 'busy');
%% extract positions and bin image
binnedImage = im2uint16(visModuleCluster(positions, 'histogramBinning', binSize));
%% extract nucleus to background-free image
% flattens intensity spikes
filteredImage = binnedImage;
for i = 1:nFilters
    h = fspecial(filterType);
    filteredImage = imfilter(filteredImage, h);
end
% binarize image by thresholding using Otsu's method
threshold = graythresh(filteredImage);
binnarizedImage = imbinarize(filteredImage, threshold);
binnarizedImage = imfill(binnarizedImage, 'holes');
% remove small objects and extract nuclei
extractedNuclei = immultiply(binnedImage, bwareaopen(binnarizedImage, maxExternalObject));
%% image reduction by a factor of PixRedFactor
inversePixRedFactor = 1/PixRedFactor;
reducedImage = imresize(extractedNuclei, inversePixRedFactor);
%% SOBEL edge detection application
sobelImage = edge(reducedImage, 'Sobel', SobelThresh);
%% remove outer edge (perimeter) of nucleus
innerNuclei = imfill(logical(reducedImage), 'holes');
% perimeter subtraction by (n)th times
for i = 1:nPerimeter
    extractedPerimeter = bwperim(innerNuclei);
    % remove perimeter i.e. subtract pixels
    innerNuclei = logical(innerNuclei-extractedPerimeter);
end
%% summarize results
% get number of connected objects and remove small residuals from perimeter
% removal
innerNuclei = bwareaopen(innerNuclei, maxExternalObjectReduced);
nObjects = bwconncomp(innerNuclei);
result = zeros(nObjects.NumObjects, 4);
areaObjects = regionprops(nObjects, 'Area');
circularityObjects = regionprops(nObjects, 'Eccentricity');
[row, col] = size(innerNuclei);
% number of edges
for ii = 1:nObjects.NumObjects
    result(ii, 1) =  areaObjects(ii).Area;
    result(ii, 2) =  circularityObjects(ii).Eccentricity;
    templateImage = zeros(row, col);
    % get n-th object
    templateImage(nObjects.PixelIdxList{1, ii}) = 1;
    % number of edges
    result(ii, 3) = sum(sum(immultiply(sobelImage, templateImage)));
    % density of edges
    result(ii, 4) = ( result(ii, 3) ./ result(ii, 1) ) * 100;
end
%% plot image
if plotOverlay == true
    visImage(sobelImage, innerNuclei, PixRedFactor, positions, binSize);
end
%% output
dataMat(1).compactionParameter = result;
if isRandom == true
    obj.randomClusterStruct = dataMat;
else
    obj.clusterStruct = dataMat;
end
multiWaitbar('CLOSEALL');
end




%% plotting function
function visImage(sobelImage, innerNuclei, PixRedFactor, positions, binSize)
    edgeImage = uint8(immultiply(sobelImage, innerNuclei))*255;
    edgeImage = imresize(edgeImage, PixRedFactor);
    SRimage = visModuleCluster(positions, 'gaussianBlur', binSize);
    figure
    subplot(1, 3, 1)
    imshow(imadjust(uint8(SRimage)))
    title('SMLM image')
    subplot(1, 3, 2)
    imshow(edgeImage)
    title('Detected Edges')
    subplot(1, 3, 3)
    imshowpair(imadjust(uint8(SRimage)), imadjust(edgeImage), 'falsecolor')
    title('Overlay')
end