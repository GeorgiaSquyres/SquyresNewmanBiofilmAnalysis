% Biofilm segmentation based on manual thresholding

% Georgia Squyres, Newman Lab, Caltech

function mask = biofilmSeg_Threshold(img)

%% Threshold

% Threshold values are specified manually, with a correction applied
% for intensity changes in Z and over time. 

threshes = linspace(175,125,sizeZ); 
tCorrection = linspace(1,1.05,sizeT);
filtimg = img;
for t = 1:sizeT
    for z = 1:sizeZ
        filtimg(:,:,z,t) = imgaussfilt(filtimg(:,:,z,t),5);
    end
end
mask = filtimg;
for t = 1:sizeT
    for z = 1:sizeZ
        currImg = filtimg(:,:,z,t);
        thresh = threshes(z).*tCorrection(t);
        m = currImg<thresh;
        m(img(:,:,z,t)==0) = 0;
        mask(:,:,z,t) = m;
    end
end

% Clean up mask: remove smaller components

mask2 = zeros(size(mask));
for t = 1:sizeT
    m = mask(:,:,:,t);
    m2 = zeros(size(m));
    cc = bwconncomp(m);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] = max(numPixels);
    m2(cc.PixelIdxList{idx}) = 1;
    mask2(:,:,:,t) = m2;
end

mask = mask2; clear mask2

% View center slice of mask over time to check segmentation
%{
for i = 1:size(mask,4)
    imagesc(squeeze(mask(:,floor(size(mask,2)/2,:,i)))
    pause(0.5)
end
%}

%% Save mask

% create image file to save
[tiffImg,tagstruct] = createTiffFile;


% save first img
for t = 1:sizeT
    for z = 1:sizeZ
        tiffImg.setTag(tagstruct)
        tiffImg.write(uint8(mask(:,:,z,t)));
        tiffImg.writeDirectory(); % saves a new page in the tiff file
    end
end

%% Child functions

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
        tiffImg = Tiff([outName,'_mask.tif'],'w8');
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

end

