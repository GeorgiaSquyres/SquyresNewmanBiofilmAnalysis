% Identifies lysis events and measures their depth in the biofilm based on 
% whole-biofilm segmentation. Biofilm segmentation can either be performed 
% using the threshold-based method provided, or masks can be loaded that 
% were generated through a separate method (e.g. Cellpose). Lysis events 
% can either be detected using automatic peak finding or, as a verification 
% method, can be manually identified in the provided interactive interface.
% To toggle between these modes, uncomment the relevant code sections.
% Lysis depths are compared to a random simulation and results are plotted.

% Georgia Squyres, Newman Lab, Caltech

function getLysisDepths

addpath('/bfmatlab/')

%% Load image
autoloadBioFormats = 1;
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

[file, path] = uigetfile(bfGetFileExtensions, 'Choose a file to open');
img = bfopen(fullfile(path,file)); % load image to segment

% parse metadata
ht = img{2};
sizeX = str2double(ht.get('Global SizeX'));
sizeY = str2double(ht.get('Global SizeY'));
sizeZ = str2double(ht.get('Global SizeZ'));
sizeT = str2double(ht.get('Global SizeT'));
MD = img{4};
pixelSizeX = double(MD.getPixelsPhysicalSizeX(0).value()); % in µm
pixelSizeY = double(MD.getPixelsPhysicalSizeY(0).value()); % in µm
pixelSizeZ = double(MD.getPixelsPhysicalSizeZ(0).value()); % in µm
pixelSizeT = 240; % imaging time interval in minutes

% concatenate to 4D image
img2 = [];
for i = 1:sizeT
    temp = cat(3,img{1}{(i-1)*sizeZ+1:i*sizeZ,1});
    img2 = cat(4,img2,temp);
end
img = img2; clear img2; clear temp;

%% Segment biofilm by thresholding

mask = biofilmSeg_Threshold(img);

%% Alternately, load mask from external segmentation method
%{
[file, path] = uigetfile(bfGetFileExtensions, 'Open mask file');
mask = bfopen(fullfile(path,file)); % load mask
% concatenate to 4D image
img2 = [];
for i = 1:sizeT
    temp = cat(3,mask{1}{(i-1)*sizeZ+1:i*sizeZ,1});
    img2 = cat(4,img2,temp);
end
mask = img2; clear img2; clear temp;
%}

%% Get biofilm surface from mask

perim = zeros(size(mask));
for t = 1:size(mask,4)
    for x = 1:size(mask,1)
        for y = 1:size(mask,2)
            currInt = mask(x,y,:,t);
            zCoord = find(currInt==1);
        end
    end
    % Remove bottom edge pixels
    perim(:,:,:,t) = bwperim(mask(:,:,:,t));
    perim(1,:,:,t) = 0;
    perim(end,:,:,t) = 0;
    perim(:,1,:,t) = 0;
    perim(:,end,:,t) = 0;
    perim(:,:,1,t) = 0;
    perim(:,:,end,t) = 0;
end

%% Find points automatically by peak finding

pointsList = lysisCaller_Automated(img);

%% Alternately, manually annotate points in GUI

%{
pointsList = lysisCaller_Interactive(img);
%}

%% Compute depths from points list

out = 0;
depths = [];
times = [];
for t = 1:sizeT
    currPerim = perim(:,:,:,t);
    ind = find(currPerim==1);
    [x,y,z] = ind2sub(size(currPerim),ind);
    currPeaks = pointsList(pointsList(:,4)==t,:);
    for p = 1:size(currPeaks,1)
        if mask(round(currPeaks(p,2)),round(currPeaks(p,1)),round(currPeaks(p,3)),t) == 1 
            dists = pdist2([x y z],currPeaks(p,[2 1 3]));
            depths(end+1) = min(dists);
            times(end+1) = t;
        else
            out = out+1;
        end
    end
end

%% simulate

allSimDepths = [];
for t = 1:sizeT
    currPerim = perim(:,:,:,t);
    ind = find(currPerim==1);
    [x1,y1,z1] = ind2sub(size(currPerim),ind);
    numPoints = sum(times==t);
    numPoints = numPoints*100;
    x = ceil(rand(numPoints*10,1).*sizeX);
    y = ceil(rand(numPoints*10,1).*sizeY);
    z = ceil(rand(numPoints*10,1).*sizeZ);
    simDepths = [];
    for i = 1:length(x)
        if mask(x(i),y(i),z(i),t) == 1 && z(i) >= 4
            simDepths(end+1) = min(pdist2([x(i) y(i) z(i)],[x1 y1 z1]));
        end
        if length(simDepths) == numPoints
            break
        end
    end
    allSimDepths = [allSimDepths,simDepths];
end

%% Plot results

[n,edges] = histcounts(allSimDepths,'Normalization','probability','BinWidth',1);
plot(mean([edges(1:end-1);edges(2:end)]),n,'k--','LineWidth',2); hold on;
histogram(depths,'Normalization','probability','BinWidth',1,'FaceColor',[0    0.4470    0.7410]);
xlim([0 22]);
set(gca,'LineWidth',2,'FontSize',18,'TickDir','out','Box','off')
axis square;

