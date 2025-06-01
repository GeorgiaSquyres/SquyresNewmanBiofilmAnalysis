% Automated identification of cell lysis events based on peak finding. 
% Note that peak intensity values are manually specified. 

% Georgia Squyres, Newman Lab, Caltech

function [peaksList] = lysisCaller_Automated(img)

%% Find and filter peaks, frame by frame

peaksList = [];

% user specified spot intensity: values should decrease w/ increasing Z
intCutStart = 220; 
intCutEnd = 20; 
intensityThreshes = linspace(intCutStart,intCutEnd,size(img,3));

fprintf('Finding peaks: ');
for i = 1:size(img,4)
    fprintf('.');
    
    % select current frame
    currImg = img(:,:,:,i);
    % replace NaNs with zeros
    currImg(isnan(currImg))=0;
    
    % 3D peak finding
    bw = imregionalmax(imgaussfilt3(currImg,2)); % gaussian filter applied
    peaks = find(bw==1);
    [x,y,z] = ind2sub(size(currImg),peaks);
    
    % measure intensity at each 3D peak
    intensities = zeros(size(x));
    for ii = 1:length(x)
        intensities(ii) = currImg(x(ii),y(ii),z(ii));
    end
    
    % intensity filter, variable with depth
    for ii = 1:length(x)
        intensityThresh = intensityThreshes(z(ii));
        if intensities(ii) >= intensityThresh
            peaksList = [peaksList;[x(ii) y(ii) z(ii) i]];
        end
    end
end

fprintf('\n');

%{
% show peaks on image
figure; 
imagesc(max(img(:,:,:,end),[],3)); hold on; plot(peaksList(:,2),peaksList(:,1),'.r');
figure; 
imagesc(squeeze(max(img(:,:,:,end),[],2))); hold on; plot(peaksList(:,3),peaksList(:,1),'.r');
%}

%% Peak linking

mergeDist = 5; 
windowSize = 1000; % for memory management
doneMerging = false; doneCounter = 0;
weightList = ones(length(peaksList),1); % weights for weighted averaging
loopCounter = 999;
while ~doneMerging
    % compute all euclidian distances
    loopCounter = loopCounter+1;
    if loopCounter == 1000
        disp(length(peaksList));
        loopCounter = 0;
    end
    % shuffle points order
    peaksList = peaksList(randperm(length(peaksList)),:);
    % currDistances(i,j) = euclidian distance between i and j
    if windowSize > length(peaksList)
        windowSize = length(peaksList);
    end
    currDistances = squareform(pdist(peaksList(1:windowSize,1:3)));
	% end if minimum distance is greater than merge cutoff
    minDistance = min(currDistances(currDistances>0));
    if size(peaksList,1) == 1 || minDistance > mergeDist
        if windowSize == length(peaksList) || doneCounter == 100
            doneMerging = true;
        else
            doneCounter = doneCounter + 1;
        end
    else
        % find indices of pair to be merged
        doneCounter = 0; 
        [i,j] = find(currDistances==minDistance);
        i = i(1); j = j(1); % fix when multiple pairs have the same distance
        % get new set of average coordinates 
        newWeight = (weightList(i)+weightList(j));
        newCoords = (peaksList(i,1:3)*weightList(i)+peaksList(j,1:3)*weightList(j))./newWeight;
        % add earlier of two time coordinates
        newCoords = [newCoords,min(peaksList(i,4),peaksList(j,4))]; 
        % remove old spots from list
        peaksList([i,j],:) = [];
        weightList([i,j],:) = [];
        % add new merged spot to list
        peaksList(end+1,:) = newCoords;
        weightList(end+1,:) = newWeight;
    end
end

% Show lysis events on time lapse
%{
figure; 
for i = 1:size(img,4)
    currFrame = img(:,:,:,i);
    imshow(squeeze(max(currFrame,[],2)),[],'InitialMagnification',200);
    hold on;
    currPeaks = peaksList(peaksList(:,4)==i,:);
    plot(currPeaks(:,3),currPeaks(:,1),'.r','MarkerSize',10);
    hold off;
    pause(0.5);
end
%}
