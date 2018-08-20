function [rois,regions]=kyoMakeROIsFromCCimage(ccimage,minArea,displayit,override)

maxArea = minArea*10;
if nargin<4
    a=[100 100];
    b=-0.1;
else
    a=override(1:2);
    b=override(3);
end

bw=adaptivethreshold(ccimage,a,b,0); % numbers determined emperically

if nargin<2
    displayit=0;
end

regions=bwlabeln(bw);
s=regionprops(regions);

roiNum=1;
for r=1:max(regions(:))
   if s(r,1).Area<minArea
       regions(regions==r)=0;
   elseif s(r,1).Area>maxArea
       regions(regions==r)=0;
   else
       rois(:,:,roiNum)=double(regions==r);
       roiNum=roiNum+1;
   end
end

regions_th=regions>0;
regions=bwlabeln(regions_th);

if displayit==1
    figure('units','normalized','outerposition',[0.4 0.6 .5 .3])
    subplot(1,3,1)
    imagesc(ccimage)
%     axis square off
    subplot(1,3,2)
    imagesc(bw)
%     axis square off
    subplot(1,3,3)
    imagesc(regions)
%     axis square off
end
