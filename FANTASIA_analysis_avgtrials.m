% Combine several trials and make a image stack with lower SNR
clear all, clc

para.date = '180625';
para.badframes = 10;
para.fs = 23.49;
para.imgsize = [653 652];
imgfolder = ['D:\=data=\80Z_imaging\img_2p\M80Z_',para.date,'\registered_NoRMCorre2'];
roipath = ['D:\=data=\80Z_imaging\img_2p\M80Z_',para.date,'_analysis'];
addpath(genpath(imgfolder))
% for tonepips session
para.stim_start = 2.7 - para.badframes / para.fs; % unit: second
para.stim_end = 14.6 + para.stim_start;
para.t_total = 20;

% for music session
% para.stim_start = 4 - para.badframes / para.fs; % unit: second
% para.stim_end = 102 + para.stim_start;
% para.t_total = 110;


%% load image stacks and extract Ca trace 
imgfolderlist = dir(fullfile(imgfolder,'*.tif'));
filerange = [1 4 5 9];
I_avg = zeros([para.imgsize, 2574]);
tic 
for i = 1:length(filerange)
    filename = imgfolderlist(filerange(i)).name
% extract Ca traces from this image stack, and save it as a matrix: 
% each column is the trace from one ROI
    I_stack = loadtiff(fullfile(imgfolder,filename));
    I_stack = double(I_stack);
    I_avg = (1/i).*I_stack + ((i-1)/i).*I_avg;
end
toc

imshow(I_avg(:,:,40),[])
I_avg = uint16(I_avg);
saveastiff(I_avg, fullfile(roipath,'avg_movie_music4trials.tif'))


%% Convert image stack to relative values
I_stack = loadtiff(fullfile(roipath,'TrialAvg_ToneDown_13-22.tif'));
I_stack = double(I_stack);
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
I_stack_rel = ((I_stack_rel - minmatCrix)./(maxmatrix - minmatrix)) .* (2^16-1);
savename = fullfile(roipath,'TrialAvg_ToneUp_2-12_relativeafter.tif');
saveastiff(uint16(I_stack_rel),savename)


