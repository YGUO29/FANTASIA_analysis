% temporary data process code
addpath(genpath(cd))
clear all, clc
listing = dir('M80Z_180625/Ca_trace_ToneUp/*.mat');
nTrial = length(listing);

Ca_trace_ToneUp = cell(nTrial,2);
for i = 1:nTrial
    load(listing(i).name)
    Ca_trace_ToneUp{i,1} = listing(i).name;
    Ca_trace_ToneUp{i,2} = Ca_trace;
end

%% temp
% i = 10;
% Ca_trace_ToneDown{i,1} = 'Ca_traces_180625_4144.mat';
% Ca_trace_ToneDown{i,2} = temp(:,2:end);
% save('D:\=data=\80Z_imaging\img_2p\M80Z_180625\Ca_trace_ToneDown_wBG','Ca_trace_ToneDown');
%% load structure (multi-trial)
addpath(genpath(cd))
clear all, clc
load Ca_trace_ToneUp_wBG.mat
% 1st column = session name, 2nd column = Ca traces (nFrame x nUnit)
nFrame = size(Ca_trace{1,2},1);
nUnit = size(Ca_trace{1,2},2);
nTrial = size(Ca_trace,1);
fStart = 63;                        % frame number before and after stim (~2.7s)
fLength = 470;                      % frame number for trial duration (~20s)
fEnd = fLength - fStart;            % frame number for stim end


%% Calculate averaged trace
trace_all = Ca_trace(:,2);
trace_all = cell2mat(trace_all');   %nFrame x (nUnit*nTrial)
trace_avg = zeros( nFrame, nUnit);
for i = 1:nUnit
    trace_avg(:,i) = mean(trace_all(:,i:nTrial:end),2);
end
%%  Calculate relative values
trace = trace_avg;                % use avaraged trace OR
% trial = 2; trace = Ca_trace{trial,2};      % use single trial trace (first dimension = trial number)
Ca_base = zeros(size(trace,2),1);

for i = 1:length(Ca_base)
    Ca_base(i) = mean(trace(10:fStart,i));
end

% deltaF/F calculation
trace_rel = trace;
for i = 1:length(Ca_base)
    trace_rel(:,i) = (trace(:,i) - Ca_base(i).*ones(size(trace,1),1))./Ca_base(i);
end

%% figures (strong unit #1 and dendrites)
% i = [2 9 11 16 19 30];
figure,
i = [1 13];
maxAmp = max(max(trace_rel(:,i)));
minAmp = min(min(trace_rel(11:end,i)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);
plot(trace_rel(:,i),'linewidth',2),title(['Trial ', num2str(trial), ', Unit ',num2str(i)]),hold on
% plot(trace_rel(:,i),'linewidth',2),title(['Average', ', Unit ',num2str(i)]),hold on

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
%% figures (all others)
% unit 5 is eliminated, since it has too low baseline and the resulting
% deltaF/F is huge
i = [3 4 5 6 7 8 9 10 11 12];
figure,
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
%% figures (whole fov)
i = 16;
figure,
maxAmp = max(max(trace_rel(:,i)));
minAmp = min(min(trace_rel(11:end,i)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
% plot(trace_rel(:,i),'linewidth',2),title(['Unit ',num2str(i)]),hold on

xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);
% plot(trace_rel(:,i),'linewidth',2),title(['Trial ', num2str(trial), ', Entire FOV']),hold on
plot(trace_rel(:,i),'linewidth',2),title(['Average, Entire FOV']),hold on
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
%% figures (plot separately)
% index = [2 9 11 16 19 30];
index = [1 2 13 14 15];
nUnit = length(index);
CM = lines(nUnit);
maxAmp = max(max(trace_rel(:,index)));
figure,
for i = 1:nUnit
subplot(nUnit,1,i),
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
plot(trace_rel(:,index(i)),'linewidth',3,'color',CM(i,:))

axis([10 fLength -0.4 maxAmp])
set(gca,'xtick',[],'ytick',[]);
set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1])
ylabel(['unit ',num2str(index(i))],'fontsize',20,'color','k'),set(get(gca,'YLabel'),'Rotation',0)
% title(['Unit ',num2str(i)])
% set(gca,'Visible','off')
% obj = scalebar;
if i == 1
    hold on,line([fEnd+20 fEnd+20],[0.7 1.7],'color','k','linewidth',3)
    hold on,line([fEnd+19 fEnd+19+23.49],[0.7 0.7],'color','k','linewidth',3)
end
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3,.8]);


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

%% Analyze background etc.
trial = 1;      % which trial
i = 1:17;       % which unit

trace_rel = Ca_trace{trial,2};
% trace_rel = trace_avg;

figure,
maxAmp = max(max(trace_rel(:,i)));
minAmp = min(min(trace_rel(11:end,i)));
axis([10 fLength minAmp-0.1 maxAmp])
hold on, area([fStart fEnd],[maxAmp maxAmp],'facecolor',0.9.*ones(1,3),'linestyle','none')
hold on, line([0 fLength],[0 0],'linewidth',2,'color',[.5 .5 .5],'linestyle','--')
% plot trace from units
plot(trace_rel(:,i(1:end-2)),'linewidth',2),title(['Unit ',num2str(i)]),hold on
% plot trace from whole FOV
plot(trace_rel(:,i(end-1)),'linewidth',3,'color','m'),title(['Unit ',num2str(i)]),hold on
% plot trace from background
plot(trace_rel(:,i(end)),'linewidth',3,'color','k'),title(['Unit ',num2str(i)]),hold on

xticks(fStart:23.49:fEnd)
xticklabels({''})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, .4, .2]);