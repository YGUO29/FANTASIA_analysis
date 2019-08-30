% extract Ca trace from imageJ rois
clear all, clc

para.date = '180419';
para.badframes = 10;
para.fs = 23.49;
para.imgsize = [653 652];
imgfolder = ['D:\=data=\80Z_imaging\img_2p\M80Z_',para.date,'\registered_NoRMCorre2'];
addpath(genpath(imgfolder))
% for tonepips session
% para.stim_start = 2.7 - para.badframes / para.fs; % unit: second
% para.stim_end = 14.6 + para.stim_start;
% para.t_total = 20;

% for music session
para.stim_start = 4 - para.badframes / para.fs; % unit: second
para.stim_end = 102 + para.stim_start;
para.t_total = 110;

%% load ROI file
% roi set file path:
roipath = ['D:\=data=\80Z_imaging\img_2p\M80Z_', para.date,'_analysis'];
roiname = [roipath,'\RoiSet2.zip'];
% read roi to MATLAB, with ReadImageJROI.m
addpath(genpath('D:\=code=\Toolbox'))
[sROI] = ReadImageJROI(roiname);

%% create ROI masks
masks = cell(1,length(sROI));
masks_outline = masks;
for i = 1:length(sROI)
    masks{i} = logical(zeros(para.imgsize));

    switch sROI{i}.strType
        case 'Freehand'
            sROI{i}.ind_linear = sub2ind(para.imgsize,sROI{i}.mnCoordinates(:,2),sROI{i}.mnCoordinates(:,1));

        case 'Rectangle'
%         sROI.vnRectBounds: ['nTop', 'nLeft', 'nBottom', 'nRight']    
            x = [sROI{i}.vnRectBounds(2),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(4),sROI{i}.vnRectBounds(2)]; % 4 points in clockwise order
            y = [sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(1),sROI{i}.vnRectBounds(3),sROI{i}.vnRectBounds(3)]; % 4 points in clockwise order
            sROI{i}.bw =poly2mask(x, y, para.imgsize(1), para.imgsize(2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            sROI{i}.ind_linear = sub2ind(para.imgsize,sROI{i}.b(:,1),sROI{i}.b(:,2));

        case 'Polygon'
            x = sROI{i}.mnCoordinates(:,1);
            y = sROI{i}.mnCoordinates(:,2);
            sROI{i}.bw =poly2mask(x, y, para.imgsize(1), para.imgsize(2)); % an area mask for polygon shape
            sROI{i}.b = bwboundaries(sROI{i}.bw); % convert areas to boundaries
            sROI{i}.b = sROI{i}.b{1};
            
            sROI{i}.ind_linear = sub2ind(para.imgsize,sROI{i}.b(:,1),sROI{i}.b(:,2));
        case 'Oval'
            
    end
    masks{i}(sROI{i}.ind_linear) = 1;
    masks{i} = imfill(masks{i},'holes');
    masks{i} = sparse(masks{i});
%     imshow(full(masks{i})), pause
end
% clear sROI
%% load image stacks and extract Ca trace 
imgfolderlist = dir(fullfile(imgfolder,'*.tif'));
filerange = 1;
traces = cell(1,length(filerange)); % each cell includes the traces from that trial
tic 
for i_trial = 1:length(filerange)
    filename = imgfolderlist(filerange(i_trial)).name
% extract Ca traces from this image stack, and save it as a matrix: 
% each column is the trace from one ROI
    I_stack = loadtiff(fullfile(imgfolder,filename));
    I_mean = mean(I_stack,3);
    % figure,imshow(I_mean,[]);

    traces{i_trial} = zeros(size(I_stack,3), length(masks)); % each column is the trace from one ROI
    for i_mask = 1:length(masks)
        for i_img = 1:size(I_stack,3)
            imgtemp = I_stack(:,:,i_img);
            traces{i_trial}(i_img,i_mask) = mean(imgtemp(masks{i_mask}));
        end
    end
end
toc

% compute deltaF/F
traces_rel = traces;
for i_trial = 1:length(filerange)

    trace = traces{i_trial};
    % ====================================
    trace_base = zeros(1, size(trace,2));
    for i = 1:length(trace_base)
        trace_base(i) = mean(trace(1:floor(para.fs*para.stim_start),i));
    end
    % deltaF/F calculation
    trace_rel = trace;
    for i = 1:length(trace_base)
        trace_rel(:,i) = (trace(:,i) - trace_base(i).*ones(size(trace,1),1))./trace_base(i);
    end
    traces_rel{i_trial} = trace_rel;
end

%% Convert image stack to relative values

para.zsize = size(I_stack,3);
I_stack_rel = zeros(size(I_stack));

% method 1: 
I_base = mean(I_stack(:,:,1:floor(para.fs*para.stim_start)),3);
I_base = repmat(I_base,[1,1,size(I_stack,3)]);
I_stack_rel = (I_stack - I_base) ./ I_base;
% method 2:,Ooo
% I_stack_rel = I_stack ./ I_base;

% method 3:
% I_base = mean(reshape(I_stack,[para.imgsize(1)*para.imgsize(2),para.zsize]));
% for i = 1:para.zsize
%     base_temp = I_base(i).*ones(para.imgsize);
%     I_stack_rel(:,:,i) = (I_stack(:,:,i) - base_temp) ./ base_temp;
% end

% convert double to uint16
maxvalue = max(I_stack_rel(:)); maxmatrix = maxvalue.*ones(size(I_stack_rel));
minvalue = min(I_stack_rel(:)); minmatrix = minvalue.*ones(size(I_stack_rel));
I_stack_rel = ((I_stack_rel - minmatrix)./(maxmatrix - minmatrix)) .* (2^16-1);
savename = fullfile(roipath,'4137_relative_stack.tif');
saveastiff(uint16(I_stack_rel),savename)
%%
save(fullfile(roipath,'5855-5959.mat'), 'traces', 'traces_rel', 'I_mean');

%% plot all trials for all ROIs

figure, 
set(gcf,'color','white')

traces_mat = cell2mat(traces_rel);
nROI = length(masks);
nTrial = length(traces_rel);
for i = 1:nROI % i: ROI index
    traces_temp = traces_mat(:,i:nROI:i+nROI*(nTrial-1));
    traces_mean = mean(traces_temp,2);
    clf
    % Ca trace plot    
    subplot(1,3,[1,2])
    t = para.badframes/para.fs : 1/para.fs : para.t_total;
    maxAmp = max(max(traces_temp));
    minAmp = min(min(traces_temp));
    axis([min(t) max(t) minAmp-0.1 maxAmp+0.1])
    hold on, area([para.stim_start para.stim_end],[maxAmp+0.1 maxAmp+0.1],'facecolor',0.9.*ones(1,3),'linestyle','none')
    hold on, line([0 max(t)],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
    xticks(para.stim_start:2.4:para.stim_end)
    xticklabels({'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'})
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.3, .7, .6])
%     title(['Trial ', num2str(i_trial), ', Unit ',num2str(i)])
    p1 = plot(t, traces_temp, 'linewidth',1.5); hold on
    p2 = plot(t, traces_mean, 'linewidth',3,'color','k'); hold on
    
    trialstr = num2str((1:nTrial)'); labels = [repmat('Trial ', nTrial, 1), trialstr];
    legend(p1,labels);
    title('Responses to tone pips (upward)','fontsize',16)
    xlabel('time (s)','fontsize',16)
    ylabel('deltaF/F','fontsize',16)

    % Mean image + ROI plot
    subplot(1,3,3),
    I_meannorm = (I_mean - min(min(I_mean)))./(max(max(I_mean)) - min(min(I_mean)));
    I_meannorm(sROI{i}.ind_linear) = 0.7;
    imshow(I_meannorm,[])

    pause
    
end

%% plot all trials for all ROIs (for sessions without identical length)
traces_size = cellfun(@size,traces_rel,'uniformoutput',0);
traces_size = cell2mat(traces_size'); 
min_length = min(traces_size(:,1));
for i = 1:length(traces_rel)
traces_rel{i} = traces_rel{i}(1:min_length,:);
end

figure, 
traces_mat = cell2mat(traces_rel);
nROI = length(masks);
nTrial = length(traces_rel);
for i = 1:nROI % i: ROI index
    traces_temp = traces_mat(:,i:nROI:i+nROI*(nTrial-1));
    traces_mean = mean(traces_temp,2);
    clf
    % Ca trace plot    
    subplot(1,2,1)
    t = para.badframes/para.fs : 1/para.fs : para.t_total;
    maxAmp = max(max(traces_temp));
    minAmp = min(min(traces_temp));
%     axis([min(t) max(t) minAmp-0.1 maxAmp+0.1])
    hold on, area([para.stim_start para.stim_end],[maxAmp+0.1 maxAmp+0.1],'facecolor',0.9.*ones(1,3),'linestyle','none')
    hold on, line([0 max(t)],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
    xticks(para.stim_start:2.4:para.stim_end)
    xticklabels({'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'})
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.3, .7, .6])
    title(['Trial ', num2str(i_trial), ', Unit ',num2str(i)])
    p1 = plot(traces_temp, 'linewidth',1.5); hold on
    p2 = plot(traces_mean, 'linewidth',3,'color','k'); hold on
    trialstr = num2str((1:nTrial)'); labels = [repmat('Trial ', nTrial, 1), trialstr];
    legend(p1,labels);

    % Mean image + ROI plot
    subplot(1,2,2),
    I_meannorm = (I_mean - min(min(I_mean)))./(max(max(I_mean)) - min(min(I_mean)));
    I_meannorm(sROI{i}.ind_linear) = 0.7;
    imshow(I_meannorm,[])

    pause
    
end

%% plot raw trace and corresponding ROI area
i_trial = 1;

figure, 
for i = 1:length(masks) % i: ROI index
clf
% Ca trace plot    
subplot(1,2,1)
t = para.badframes/para.fs : 1/para.fs : para.t_total;
maxAmp = max(max(traces_rel{i_trial}(:,i)));
minAmp = min(min(traces_rel{i_trial}(:,i)));
axis([min(t) max(t) minAmp-0.1 maxAmp+0.1])
hold on, area([para.stim_start para.stim_end],[maxAmp+0.1 maxAmp+0.1],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 max(t)],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(para.stim_start:2.4:para.stim_end)
xticklabels({'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.3, .7, .6]);
plot(t,traces_rel{i_trial}(:,i), 'linewidth',2), title(['Trial ', num2str(i_trial), ', Unit ',num2str(i)]),hold on

% Mean image + ROI plot
subplot(1,2,2),
I_meannorm = (I_mean - min(min(I_mean)))./(max(max(I_mean)) - min(min(I_mean)));
I_meannorm(sROI{i}.ind_linear) = 0.7;
imshow(I_meannorm,[])

pause
end

%% plot all mean traces and all ROIs

figure, 
set(gcf,'color','white')
traces_mat = cell2mat(traces_rel);
nROI = length(masks);
nTrial = length(traces_rel);
traces_mean = zeros(size(traces_mat,1), nTrial);
I_meannorm = (I_mean - min(min(I_mean)))./(max(max(I_mean)) - min(min(I_mean)));
for i = 1:nROI % i: ROI index
    traces_temp = traces_mat(:,i:nROI:i+nROI*(nTrial-1));
    traces_mean(:,i) = mean(traces_temp,2);
%     if max(traces_mean(:,i))>3
    I_meannorm(sROI{i}.ind_linear) = 0.7;
%     end
end
% Ca trace plot    
subplot(1,3,[1 2])
t = para.badframes/para.fs : 1/para.fs : para.t_total;
maxAmp = max(max(traces_mean));
minAmp = min(min(traces_mean));
axis([min(t) max(t) minAmp-0.1 maxAmp+0.1])
hold on, area([para.stim_start para.stim_end],[maxAmp+0.1 maxAmp+0.1],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 max(t)],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(para.stim_start:2.4:para.stim_end)
xticklabels({'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'})
% xticklabels(flip({'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10'}))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.3, .7, .6])


p2 = plot(t,traces_mean, 'linewidth',1);
% p2 = plot(t,traces_mean(:,max(traces_mean)>3), 'linewidth',1);

ROIstr = num2str((1:nROI)'); labels = [repmat('ROI ', nROI, 1), ROIstr];
% legend(p2,labels);
title('Responses to tone pips (upward)','fontsize',16)
xlabel('time (s)','fontsize',16)
ylabel('deltaF/F','fontsize',16)

% Mean image + ROI plot
subplot(1,3,3),
imshow(I_meannorm,[])

% =========================== above: modified 9/27/2018 =================================


%% select units to plot together 
i = [2 9 11 16 19 30];
figure,
% i = [1 13];
maxAmp = max(max(traces_mat(:,i)));
minAmp = min(min(traces_mat(:,i)));
axis([0 para.t_total minAmp-0.1 maxAmp])
hold on, area([para.stim_start para.stim_end],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 para.t_total],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(para.stim_start:23.49:para.stim_end)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);
plot(t,traces_mat(:,i),'linewidth',2),title(['Trial ', num2str(trial), ', Unit ',num2str(i)]),hold on
% plot(trace_rel(:,i),'linewidth',2),title(['Average', ', Unit ',num2str(i)]),hold on

%% select units to plot separately
index = [2 9 11 16 19 30];
nUnit = length(index);
CM = lines(nUnit);
maxAmp = max(max(traces_mat(:,index)));
figure,
for i = 1:nUnit
subplot(nUnit,1,i),
hold on, area([para.stim_start para.stim_end],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 para.t_total],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
plot(t, traces_mat(:,index(i)),'linewidth',3,'color',CM(i,:))

axis([0 para.t_total -0.4 maxAmp])
set(gca,'xtick',[],'ytick',[]);
set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1])
ylabel(num2str(i),'fontsize',20,'color','k'),set(get(gca,'YLabel'),'Rotation',0)
% title(['Unit ',num2str(i)])
% set(gca,'Visible','off')
% obj = scalebar;
if i == 1 % plot the scalebar on first figure
    hold on,line([para.stim_end+1 para.stim_end+1],[1 2],'color','k','linewidth',3)
    hold on,line([para.stim_end+0.9 para.stim_end+5.9],[1 1],'color','k','linewidth',3)
end
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3,.8]);


%% figures (strong unit #2 and dendrites)
figure,
i = [2 14 15];
maxAmp = max(max(trace_rel(:,i)));
minAmp = min(min(trace_rel(11:end,i)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
% plot(trace_rel(:,i),'linewidth',2),title(['Unit ',num2str(i)]),hold on

xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);
plot(trace_rel(:,i),'linewidth',2),title(['Trial ', num2str(trial), ', Unit ',num2str(i)]),hold on
% plot(trace_rel(:,i),'linewidth',2),title(['Average', ', Unit ',num2str(i)]),hold on


%% figure (correlation matrix)
figure,
V = trace_rel(:,[1 13 2 14 15 16]);
M = corrcoef(V);
n = size(V,2);
L = {'S1','D1','S2','D2.1','D2.2','whole'};
imagesc(M); % plot the matrix
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', L); % set x-axis labels
set(gca, 'YTickLabel', L); % set y-axis labels
title(['Correlation Matrix, Trial ',num2str(trial)], 'FontSize', 14); % set title
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar



%% low-pass filter??
% Konnerth: low-pass filtering with impulse response 50ms and cutoff
% frequency was sampling frequency devided by 10
fs = 23.49;
d1 = designfilt('lowpassiir','FilterOrder',10, ...
    'HalfPowerFrequency',fs/5,'SampleRate',fs,'DesignMethod','butter');
fvtool(d1)
x = trace_rel(10:end,i);
y = filtfilt(d1,x);
figure,subplot(2,1,1),plot(x)
subplot(2,1,2),plot(y)
