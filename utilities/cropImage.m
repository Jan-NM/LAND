function [] = cropImage()
%cropImage crops SMLM image
%
% partially based on cropOrte.m by Martin Hagmann
% crops several ROIs from batch of images
%
% Jan Neumann, 06.03.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixelSize = 10; % default SR pixel size for histogram binning
[FileName, PathName] = uigetfile('*.mat', 'Select the localization data file(s).', 'MultiSelect', 'on');
if isequal(FileName, 0)
    errorMessage = sprintf('Error no file selected');
    uiwait(warndlg(errorMessage));
    return;
end

FileName = cellstr(FileName); % convert to cell if only one file is selected

for ii = 1:numel(FileName)
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1, 2, 1);
    % load list of localizations
    locFile = load(fullfile(PathName, FileName{ii}));
    oldField = char(fieldnames(locFile));
    if ~(strcmp(oldField, 'Orte'))
        newField = 'Orte';
        [locFile.(newField)] = locFile.(oldField);
        locFile = rmfield(locFile, oldField);   
    end
    % reconstruct histogram binned image for cropping
    x1 = round(locFile.Orte(:, 2)./pixelSize);
    y1 = round(locFile.Orte(:, 3)./pixelSize);
    SRimage = sparse(x1-min(x1)+1, y1-min(y1)+1, 1);
    % transform to 8 bit
    SRimage = full(SRimage);
    SRimage = SRimage > 0;
    imagesc(SRimage);
    axis image;
    k = waitforbuttonpress;
    point1 = get(gca, 'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca, 'CurrentPoint');    % button up detected
    point1 = point1(1, 1:2);              % extract x and y
    point2 = point2(1, 1:2);
    p1 = min(point1, point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    % Find the coordinates of the box.
    xCoords = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    yCoords = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    x1Cropped = round(xCoords(1));
    x2Cropped = round(xCoords(2));
    y1Cropped = round(yCoords(5));
    y2Cropped = round(yCoords(3));
    hold on
    axis manual
    plot(xCoords, yCoords, 'r-'); % redraw in dataspace units
    % Display the cropped image.
%     figure
%     croppedImage = SRimage(y1Cropped:y2Cropped, x1Cropped:x2Cropped);
%     imshow(croppedImage);
%     axis on;
%     title('Region that you defined', 'FontSize', 30);
    % transform to Orte coordinates
    minX = x1Cropped * pixelSize;
    minY = y1Cropped * pixelSize;
    maxX = x2Cropped * pixelSize;
    maxY = y2Cropped * pixelSize;
    % work first with min / max and floor / ceil then multiply with pixelsize   

    % rectangle('Position',[imRect(1), imRect(2), imRect(3), imRect(4)],...
    %    'Curvature',[0.0, 0.0], 'EdgeColor', 'r', 'LineWidth', 3,...
    %    'LineStyle','-');
    OrteCropped = locFile.Orte(locFile.Orte(:, 3) > minX & locFile.Orte(:, 3) < maxX, :);
    OrteCropped = OrteCropped(OrteCropped(:, 2) > minY & OrteCropped(:, 2) < maxY, :);
    % shift list of localization
    mincroppedX = min(OrteCropped(:, 3));
    mincroppedY = min(OrteCropped(:, 2));
    OrteCropped(:, 3) = OrteCropped(:, 3) - mincroppedX + 1;
    OrteCropped(:, 2) = OrteCropped(:, 2) - mincroppedY + 1;
    x1 = round(OrteCropped(:, 2)./pixelSize);
    y1 = round(OrteCropped(:, 3)./pixelSize);
    SRimage = sparse(x1-min(x1)+1, y1-min(y1)+1, 1);
    % transform to 8 bit
    SRimage = full(SRimage);
    SRimage = SRimage > 0;
    subplot(1, 2, 2);
    imagesc(SRimage);
    axis image;
    Orte = OrteCropped;
    save([PathName FileName{ii}(1:end-4) '_cropped.mat'], 'Orte');
end
end