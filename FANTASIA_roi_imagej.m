% analyze with ROI imported from ImageJ
% roi set file path:
roipath = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625';
roiname = [roipath,'\RoiSet.zip'];

% read roi to MATLAB, with ReadImageJROI.m
addpath(genpath('\Toolbox'))
[sROI] = ReadImageJROI(roiname);

%% import mean image
filename = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625\registered_fiji\AVG_20180625T112517_registered.tif';
I = loadtiff(filename);
figure,imshow(I,[]);

%% plot all ROI contours (grayscale)
withMean = 1;
I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));
mask = zeros(size(I));
for i = 1:length(sROI)
    switch sROI{i}.strType
        case 'Freehand'
            ind_linear = sub2ind(size(I),sROI{i}.mnCoordinates(:,2),sROI{i}.mnCoordinates(:,1));
            
        case 'Rectangle'
%         sROI.vnRectBounds: ['nTop', 'nLeft', 'nBottom', 'nRight']    
            x = [sROI{i}.vnRectBounds(2),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(2)]; % 4 points in clockwise order
            y = [sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(3),sROI{i}.vnRectBounds(3)]; % 4 points in clockwise order
            sROI{i}.bw =poly2mask(x,y,size(I,1),size(I,2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            ind_linear = sub2ind(size(I),sROI{i}.b(:,1),sROI{i}.b(:,2));
            
        case 'Polygon'
            x = sROI{i}.mnCoordinates(:,1);
            y = sROI{i}.mnCoordinates(:,2);
            sROI{i}.bw =poly2mask(x,y,size(I,1),size(I,2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            ind_linear = sub2ind(size(I),sROI{i}.b(:,1),sROI{i}.b(:,2));
            

        case 'Oval'
            
    end
    if ~withMean
        mask(ind_linear) = 1;
    else
        I_norm(ind_linear) = 0.7;
    end
    
end

figure,
if ~withMean
    imshow(mask,[]) 
else
    imshow(I_norm,[]) 
end


%% plot ROIs on mean image (colored ROIs with grey mean image)
% prepare for comining with ROIs
I_norm = (I - min(min(I)))./(max(max(I)) - min(min(I)));
img = repmat(I_norm,1,1,3); % three layers, representing hue/saturation/value in hsv
img = rgb2hsv(img);
img = reshape(img,[],3); % brightness represents the mean image
% imagesc(hsv2rgb(reshape(img,size(I,1),size(I,2),3)))
% colormap('gray');
% hold all;


for i = 1:length(sROI)
    switch sROI{i}.strType
        case 'Freehand'
            ind_linear = sub2ind(size(I),sROI{i}.mnCoordinates(:,2),sROI{i}.mnCoordinates(:,1)); % x and y are inversed in ROI extracted by ImageJ
            
        case 'Rectangle'
%         sROI.vnRectBounds: ['nTop', 'nLeft', 'nBottom', 'nRight']    
            x = [sROI{i}.vnRectBounds(2),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(2)]; % 4 points in clockwise order
            y = [sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(3),sROI{i}.vnRectBounds(3)]; % 4 points in clockwise order
            sROI{i}.bw =poly2mask(x,y,size(I,1),size(I,2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            ind_linear = sub2ind(size(I),sROI{i}.b(:,1),sROI{i}.b(:,2));
            mask(ind_linear) = 1;
        case 'Polygon'
            x = sROI{i}.mnCoordinates(:,1);
            y = sROI{i}.mnCoordinates(:,2);
            sROI{i}.bw =poly2mask(x,y,size(I,1),size(I,2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            ind_linear = sub2ind(size(I),sROI{i}.b(:,1),sROI{i}.b(:,2));
            mask(ind_linear) = 1;

        case 'Oval'
        
    end
    img(ind_linear, 1) = rand;
    img(ind_linear, 2) = 1;
    img(ind_linear, 3) = 1; 
    
end


imagesc(hsv2rgb(reshape(img,size(I,1),size(I,2),3)))
set(gcf,'color','w');
axis off;
axis square;



        
    
    
    