% Image file pre-processing for confocal Z stack time lapses
% Configured for .ims files from an Andor Dragonfly confocal
% Merges multiple files with different names and potentially different Z 
%    stack heights, combines, downsamples to 8-bit, and registers
% Memory efficient: only loads pairs of images at a time at the cost of
%    cropping registered images

% Georgia Squyres, Newman Lab, Caltech 2024
% Some code is modified from the MATLAB Bioformats toolbox:
% https://www.openmicroscopy.org/bio-formats/downloads/ 

function imagePreprocessing

%% INITIALIZE

% You will need to download these files and add them to your MATLAB path
% See dependencies in README 
addpath('/bfmatlab');
addpath('/natsortfiles');

% === USER PARAMETERS ===

sizeC = 2; % number of channels
regChannel = 2; % reference channel to use for registration
downsizeFactor = 2; %1 = full size, 2 = 50%, etc. 

path = '/path-to-image-file/'; % directory containing image files
fileHeaders = {'file_name_header_'}; % image file name header
filePositions = 0:7; % which positions to process
pixelSizeT = 240; % imaging time interval in minutes
outPath = path; % path to save outputs

% =======================

% Load bioformats
autoloadBioFormats = 1;
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

% Initialize logging
bfInitLogging();

% Silence compression type warning
warning('off','imageio:tiffutils:libtiffWarning')

tic

% Loop through positions
for position = 1:length(filePositions)
  
    currPosition = filePositions(position);

    % Update progress bar
    disp(['Processing position ',num2str(currPosition)]);

    % Identify files to open
    files = {};
    % find and sort files for current position, in order of headers
    for i = 1:length(fileHeaders)
        currHeader = fileHeaders{i};
        currHeader = [currHeader,'*F',num2str(currPosition),'*.ims'];
        f = dir([path,currHeader]);
        f = {f.name};
        f = natsortfiles(f);
        % check for zero time point at end, sort to front:
        if ~isempty(strfind(f{end},'Zyla_F'))
            temp = f(end);
            f(2:end) = f(1:end-1);
            f(1) = temp;
        end
        % append to files list
        files = [files,f];
    end
    
    % Create name to save outputs
    outName = [outPath,fileHeaders{1},'_F',num2str(currPosition)];

    % Get total number of time points- needed to create tiff save
    sizeT = 0;
    for f = 1:length(files)
        currFile = [path,files{end}];
        r = bfGetReader(currFile,0);
        sizeT = sizeT + r.getSizeT();
    end

    regMatrix = zeros(1,3); % matrix for registration offsets

    % Order to process channels: registration reference channel first
    channels = 1:sizeC;
    channels(regChannel) = [];
    channels = [regChannel,channels];

    % Loop through channels
    for c = channels
        % Update progress bar
        disp(['Channel ',num2str(c),'...'])
    
        % Update progress bar
        nbytes = fprintf(['Time point 1 of ',num2str(sizeT)]);
    
        % Get metadata and parameters for 8-bit rescaling
        % Also loads first image as firstImg
        [sizeX, sizeY, sizeZ, pixelSizeX, pixelSizeY, pixelSizeZ,...
            minSizeZ, imgMax, firstImg] = getMetadata;
    
        % assign starting image as prevImg
        prevImg = zeros(sizeX,sizeY,sizeZ,'uint8');
        prevImg(:,:,1:size(firstImg,3)) = (firstImg./imgMax).*255;
        clear firstImg
    
        % create image file to save
        [tiffImg,tagstruct] = createTiffFile;
    
        % save first img
        for z = 1:sizeZ
            tiffImg.setTag(tagstruct)
            tiffImg.write(uint8(prevImg(:,:,z)));
            tiffImg.writeDirectory(); % saves a new page in the tiff file
        end
        
        tProgress = 2; % time counter for progress bar and regmat indexing
    
        % Main loop through files
        for f = 1:length(files)
    
            currFile = [path,files{f}];
            r = bfGetReader(currFile,0);
            r.setSeries(downsizeFactor-1);
            currSizeZ = r.getSizeZ();
            currSizeT = r.getSizeT();
    
            % Skip first time point in first file: has already been processed
            if f == 1 
                timePointsList = 2:currSizeT;
            else
                timePointsList = 1:currSizeT;
            end
    
            for t = timePointsList
                % Update progress bar
                fprintf(repmat('\b',1,nbytes));
                nbytes = fprintf(['Time point ',num2str(tProgress),' of ',num2str(sizeT)]);
    
                % Load image
                img = zeros(sizeX,sizeY,sizeZ,'uint8');
                for z = 1:currSizeZ
                    plane = z + currSizeZ*(c-1)+currSizeZ*sizeC*(currSizeT-1); % which plane to read
                    currPlane = double(bfGetPlane(r, plane)); % read plane
                    % rescale image to 8-bit
                    currPlane = currPlane-90;
                    currPlane = (currPlane./imgMax).*255;
                    % add to 4D image
                    img(:,:,z) = uint8(currPlane);
                end
                
                % Register to previous image, using only the common Z range
                if c == regChannel
                    % Compute registration matrix
                    regOffsets = regmat3D(cat(4,prevImg(:,:,1:minSizeZ),img(:,:,1:minSizeZ)));
                    % Make offsets cumulative
                    regOffsets = regOffsets(2,:)+regMatrix(end,:);
                    % Append offsets to registration matrix for logging
                    regMatrix = [regMatrix;regOffsets];
                else % Otherwise use previous registration
                    regOffsets = regMatrix(tProgress,:);
                end
                % Register
                regImg = doRegistration(img,regOffsets);

                % Save
                for z = 1:sizeZ
                    tiffImg.setTag(tagstruct)
                    tiffImg.write(uint8(regImg(:,:,z)));
                    tiffImg.writeDirectory(); % saves a new page in the tiff file
                end

                prevImg = img; clear img regImg;
                tProgress = tProgress + 1;

            end

            r.close();
            
        end

        close(tiffImg);
        save([outName,'_channel',num2str(c),'_imgMax.mat'],'imgMax')
        
        fprintf('\n')

    end

    save([outName,'_regMatrix.mat'],'regMatrix')

end
    
% turn warnings back on
warning('on','all')
toc

%% Child functions

function [sizeX, sizeY, sizeZ, pixelSizeX, pixelSizeY, pixelSizeZ,...
    minSizeZ, imgMax, firstImg] = getMetadata
% Function to get image metadata, including 8-bit conversion maximum

    % First, load and read all metadata from last time point in last file

    % Get the reader for current file
    currFile = [path,files{end}];
    r = bfGetReader(currFile,0);
    r.setSeries(downsizeFactor-1);

    % Get pixel sizes  in um
    MD = r.getMetadataStore(); 
    pixelSizeX = double(MD.getPixelsPhysicalSizeX(0).value()); % in µm
    pixelSizeY = double(MD.getPixelsPhysicalSizeY(0).value()); % in µm
    pixelSizeZ = double(MD.getPixelsPhysicalSizeZ(0).value()); % in µm
    
    % Get XYZ sizes
    sizeX = r.getSizeX();
    sizeY = r.getSizeY();
    sizeZ = r.getSizeZ();
    sizeC = r.getSizeC();
    currSizeT = r.getSizeT();

    % read image for reference value for 8-bit downsampling
    lastImg = zeros(sizeY,sizeX,sizeZ);
    for zz = 1:sizeZ
        % read image
        plane = zz + sizeZ*(c-1) + sizeZ*sizeC*(currSizeT-1); % which plane to read
        lastImg(:,:,zz) = bfGetPlane(r, plane); % read plane
    end
    % rescale image to 8-bit
    lastImg = double(lastImg);
    lastImg = lastImg-90; % minimum camera value from dragonfly
    lastMax = prctile(lastImg(:),99.9); 

    r.close();

    % Second, load and read Z size from first image in first file

    % Get the reader for current file
    currFile = [path,files{1}];
    r = bfGetReader(currFile,0);
    r.setSeries(downsizeFactor-1);

    minSizeZ = r.getSizeZ();

    firstImg = zeros(sizeY,sizeX,minSizeZ);
    for zz = 1:minSizeZ
        % read image
        plane = zz + minSizeZ*(c-1); % which plane to read
        firstImg(:,:,zz) = bfGetPlane(r, plane); % read plane
    end
    % rescale image to 8-bit
    firstImg = double(firstImg);
    firstImg = firstImg-90; % minimum camera value from dragonfly
    firstMax = prctile(firstImg(:),99.9); 

    imgMax = max(firstMax,lastMax).*1.2; % maximum to rescale to, will be applied to all time points
    if imgMax < 255
        imgMax = 255;
    end

    r.close();
end

function output = regmat3D(img)
% function to obtain cumulative registration matrix for 3D images
% assumes channel order XYZT
% columns 1-3 are cumulative shifts of x, y, z respectively
% registration performed by dftreg3D

    output = zeros(size(img,4),3);
    for i = 2:size(img,4)
        % register
        input = double(img(:,:,:,i));
        prev = double(img(:,:,:,i-1));
        output(i,:) = dftreg3D(fftn(prev), fftn(input));
    end
    output = cumsum(output);
end

function [tiffImg,tagstruct] = createTiffFile
% Function to initialize Tiff file

        % Create metadata
        fiji_descr = ['ImageJ=1.52p' newline ...
                'images=' num2str(sizeZ*...
                                  sizeT) newline... 
                'channels=1' newline...
                'slices=' num2str(sizeZ) newline...
                'frames=' num2str(sizeT) newline... 
                'hyperstack=true' newline...
                'mode=grayscale' newline...  
                'loop=false' newline...  
                'min=0.0' newline...      
                'max=255.0' newline...
                'spacing=' num2str(pixelSizeZ) newline...
                'finterval=' num2str(pixelSizeT) newline...
                'unit=micron'];
        
        % Create tiff writer
        tiffImg = Tiff([outName,'_channel',num2str(c),'.tif'],'w8');
        tagstruct.ImageLength = sizeX;
        tagstruct.ImageWidth = sizeY;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 8;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.Deflate;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.XResolution = 1/(pixelSizeX*downsizeFactor);
        tagstruct.YResolution = 1/(pixelSizeY*downsizeFactor);
        tagstruct.ImageDescription = fiji_descr;
end

function img = doRegistration(img,regOffsets)

    % Dimension 1
    sh = regOffsets(1);
    size1 = size(img,1); 
    if sh > 0
        img = padarray(img, [sh 0 0], 0, 'pre');
        img = img(1:size1,:,:);
    elseif sh < 0
        img = padarray(img, [-sh 0 0], 0, 'post');
        img = img(end-size1+1:end,:,:);
    end

    % Dimension 2
    sh = regOffsets(2);
    size2 = size(img,2); 
    if sh > 0
        img = padarray(img, [0 sh 0], 0, 'pre');
        img = img(:,1:size2,:);
    elseif sh < 0
        img = padarray(img, [0 -sh 0], 0, 'post');
        img = img(:,end-size2+1:end,:);
    end

    % Dimension 3
    sh = regOffsets(3);
    size3 = size(img,3); 
    if sh > 0
        img = padarray(img, [0 0 sh], 0, 'pre');
        img = img(:,:,1:size3);
    elseif sh < 0
        img = padarray(img, [0 0 -sh], 0, 'post');
        img = img(:,:,end-size3+1:end);
    end

end

end