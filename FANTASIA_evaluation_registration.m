% process binary file from suite2p, and save registered tiffs
% load regops file manually
iplane = 1; % which plane of the recording
ops = ops1{iplane};
fid  = fopen(ops.RegFile, 'r'); % opens the registered binary
[Ly Lx] = size(ops.mimg1);    % size of binary in x,y

clf;
nt0 = 0;
nimgbatch = 3000; % how many frames you can load at a time (will depend on RAM)
while 1
  data = fread(fid, Ly*Lx*nimgbatch, '*int16');
  if isempty(data)
    break;
  end
  data = reshape(data, Ly, Lx, []);
  data = data(ops.yrange, ops.xrange, :);
  NT   = size(data,3);
%   for j = 1:NT
%     imagesc(data(:,:,j),[0 2000]);%,[1000 3000]);
%     title(sprintf('frame %d', j+nt0));
%     axis square;
%     drawnow;
%     
%     pause(0.25);
%   end
%   nt0 = nt0 + NT;
end

saveastiff(data,'D:\=data=\80Z_imaging\img_2p\Suite2P_Results\80Z_imaging\img_2p\2\20180503T121108_nonrigid.tif');
%%
% evaluate motion correction results
clear all,clc
addpath('D:\=code=\NoRMCorre')

% ================= load Suite2p reesults ===========
% img2 = loadtiff('D:\=data=\80Z_imaging\img_2p\Suite2P_Results\80Z_imaging\img_2p\2\20180503T121108_nonrigid.tif');
load('D:\=data=\80Z_imaging\img_2p\Suite2P_Results\80Z_imaging\img_2p\2\regops_80Z_imaging_img_2p');
ops = ops1{1};
% ===================================================

% ================= load NoRMCorre reesults ===========
img2 = loadtiff('D:\=data=\80Z_imaging\img_2p\CaImAn_Results\20180503T121108_nonrigid_nonbidir.tif');
img2 = img2(ops.yrange,ops.xrange,:);
% ===================================================

% ================= load TurboReg reesults ===========
% img2 = loadtiff('D:\=data=\80Z_imaging\img_2p\M80Z_180419\20180419T164137_registered.tif');
% img2 = img2(ops.yrange,ops.xrange,:);
% ===================================================

% load reference image
mimg = loadtiff('D:\=data=\80Z_imaging\img_2p\M80Z_180503\AVG_20180503T121108_registered.tif');
mimg = mimg(ops.yrange,ops.xrange,:);
% method 2 - crispness calculation
mimg = mean(mimg,3);
[mimg_g,~] = imgradient(mimg);
crispness_mimg = norm(mimg_g);
%% method 1 - correlation with reference image
% load image 2 (registered image)

nFrame = size(img2,3);
C2 = zeros(nFrame,1);
for i = 1:nFrame
    C2(i) = corr2(img2(:,:,i), mimg);
end
figure,plot(C2),hold on

% method 2 - crispness calculation
img2 = mean(img2,3);
[img2_g,~] = imgradient(img2);
crispness_img2 = norm(img2_g);
%% plot correlation of unregistered image with template
clear img2
% load image 1 (unregistered image)
img1 = loadtiff('D:\=data=\80Z_imaging\img_2p\M80Z_180503\20180503T121108.tif');
img1 = img1(ops.yrange,ops.xrange,:);
nFrame = size(img1,3);
C1 = zeros(nFrame,1);
for i = 1:nFrame
    C1(i) = corr2(img1(:,:,i), mimg);
end
hold on, plot(C1)

% method 2 - crispness calculation
img1 = mean(img1,3);
[img1_g,~] = imgradient(img1);
crispness_img1 = norm(img1_g);