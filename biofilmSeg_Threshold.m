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

end

