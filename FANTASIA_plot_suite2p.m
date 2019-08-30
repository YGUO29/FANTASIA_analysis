% FANTASIA analysis - suite2p results
%% plot ROIs with mean image 

clear all, clc
filepath = 'D:\=data=\80Z_imaging\img_2p\Suite2P_Results\80Z_imaging\img_2p\3\';
filename = 'F_80Z_imaging_img_2p_plane1';
dat = load([filepath,filename]);

mimg = dat.ops.mimg1;
stat = [];
xL = [];
yL = [];
% ============ get mean image ============
mimg = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);

mimg = mimg - min(mimg(:));
mimg = mimg / max(mimg(:));

stat = cat(2, stat, dat.stat);
yL = numel(dat.ops.yrange);
xL = numel(dat.ops.xrange);

img = repmat(mimg*2,1,1,3);
img = rgb2hsv(img);
img = reshape(img,[],3);
colormap('gray');
hold all;
% ============ extract roi borders ============
% test = zeros(size(mimg)); 
for ij = 1:length(dat.stat)
%     if stat(ij).iscell
    % find pixels that are exterior to the cell 
    idist  = sqrt(bsxfun(@minus, stat(ij).xpix', stat(ij).xpix).^2 + ...
        bsxfun(@minus, stat(ij).ypix', stat(ij).ypix).^2);

    idist  = idist - diag(NaN*diag(idist));
    extpix = sum(idist <= sqrt(2)) <= 6;
    xext = stat(ij).xpix(extpix);
    yext = stat(ij).ypix(extpix);

    % color the exterior pixels
    ipix = sub2ind([yL xL], yext, xext);
%     test(ipix) = 1;
    img(ipix, 1) = rand;
    img(ipix, 2) = 1;
    img(ipix, 3) = 1; 
%     end
end

% ============ plot of one plane with circled ROIs ============
clf;
set(gcf,'color','w');
iplane = 1;
imgj = img;
imagesc(hsv2rgb(reshape(imgj,yL(iplane),xL(iplane),3)))
axis off;
axis square;

%% plot Ca traces

clear all, close all, clc
filepath = 'D:\=data=\80Z_imaging\img_2p\Suite2P_Results\80Z_imaging\img_2p\1\';
filename = 'F_80Z_imaging_img_2p_plane1_proc';
load([filepath,filename]);

traceCell = []; 
tracePil = [];
traceSpk = [];
i = 1;
for k = 1:size(dat.Fcell{1})
    if dat.stat(k).iscell
    traceCell(i,:) = dat.Fcell{1}(k,:);
    tracePil(i,:) = dat.FcellNeu{1}(k,:);
    traceSpk(i,:) = dat.sp{1}(k,:);
    i = i + 1;
    end
end

nCell = size(traceCell,1);
neuropilCoefMat = zeros(nCell,1);
for k = 1:nCell
   neuropilCoefMat(k) = dat.stat(k).neuropilCoefficient; 
end
neuropilCoefMat = repmat(neuropilCoefMat,1,size(traceCell,2));
traceRel = traceCell - neuropilCoefMat.*tracePil;


% ===== show traces =====
figure,
subplot(4,1,2),imagesc(traceCell),colorbar,title('Raw Ca trace','fontsize',12)
subplot(4,1,1),imagesc(tracePil),colorbar, title('Neuropil trace','fontsize',12)
subplot(4,1,3),imagesc(traceRel),colorbar, title('Neuropil-subtraction trace','fontsize',12)
subplot(4,1,4),imagesc(traceSpk),colorbar, title('Estimated spike rate','fontsize',12)
colormap('gray')
% ======================
fStart = 63;                        % frame number before and after stim (~2.7s)
fLength = 470;                      % frame number for trial duration (~20s)
fEnd = fLength - fStart;            % frame number for stim end

mimg = dat.ops.mimg1;
stat = [];
xL = [];
yL = [];
% ============ get mean image ============
mimg = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);

mimg = mimg - min(mimg(:));
mimg = mimg / max(mimg(:));

stat = cat(2, stat, dat.stat);
yL = numel(dat.ops.yrange);
xL = numel(dat.ops.xrange);
% colormap('gray');
% hold all;

%% plot Ca trace + Neuropil trace, Corrected trace, cell position
close all
for i = 7 % number of cell
figure,
maxAmp = max(max(traceCell(i,:)));
minAmp = min(min(tracePil(i,11:end)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);
plot(traceCell(i,:),'linewidth',2),title(['Unit ',num2str(i),', IsCell = ',num2str(dat.stat(i).iscell)]),hold on
plot(tracePil(i,:),'linewidth',2,'color','r')

% plot relative fluorescence values
figure,
maxAmp = max(max(traceRel(i,:)));
minAmp = min(min(traceRel(i,11:end)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.7, .4, .2]);
plot(traceRel(i,:),'linewidth',2,'color','k'),title(['Unit ',num2str(i),', IsCell = ',num2str(dat.stat(i).iscell)]),hold on

% plot estimated spike rates
figure,
maxAmp = max(max(traceSpk(i,:)));
minAmp = min(min(traceSpk(i,11:end)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.3, .4, .2]);
plot(traceSpk(i,:),'linewidth',2,'color','m'),title(['Unit ',num2str(i),', IsCell = ',num2str(dat.stat(i).iscell)]),hold on


% ============ extract roi borders ============
figure,
img = repmat(mimg*2,1,1,3);
img = rgb2hsv(img);
img = reshape(img,[],3);
% find pixels that are exterior to the cell 
idist  = sqrt(bsxfun(@minus, stat(i).xpix', stat(i).xpix).^2 + ...
    bsxfun(@minus, stat(i).ypix', stat(i).ypix).^2);
idist  = idist - diag(NaN*diag(idist));
extpix = sum(idist <= sqrt(2)) <= 6;
xext = stat(i).xpix(extpix);
yext = stat(i).ypix(extpix);
% color the exterior pixels
ipix = sub2ind([yL xL], yext, xext);
img(ipix, 1) = .3;
img(ipix, 2) = 1;
img(ipix, 3) = 1; 
% ============ plot of one plane with circled ROIs ============
clf; set(gcf,'color','w');
iplane = 1;

imgj = img;
imagesc(hsv2rgb(reshape(imgj,yL(iplane),xL(iplane),3)))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [.4, .3, .3, .6]);
axis off;
axis square;
end