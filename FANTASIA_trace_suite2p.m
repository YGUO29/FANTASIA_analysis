% Ca trace extraction from registered images
% by Yueqi Guo 2018/8
% based on code from Spencer Smith lab (by Yiyi Yu)
% sync
addpath(genpath('D:\=code=\Suite2P'))
addpath('D:\=code=\CaImgProc\Kyo')
datafolder = 'D:\=data=\80Z_imaging\img_2p\M80Z_18c0625';
imgfolder = 'D:\=data=\80Z_imaging\img_2p\M80Z_180625\registered';

%%  ROI by kurtosis map
folder = fullfile(imgfolder); 
SegmentPara = [20 20 -0.015]; % adaptive thresholding parameters
subtractPara = [2,5]; % neuropil correction parameters
minArea = 20; % minimal cell size (was 50)
displayit = 1; % display results
G = fspecial('gaussian',[4 4],2); % gaussian filter parameters
subsample = 5; % subsampling for gaussian filtering
imgfolderlist = dir(fullfile(imgfolder,'*_reg.bin'));
for i = 1:length(imgfolderlist)
    i
    % load image
    filename = strtrim(imgfolderlist(i).name(1:end-8))
    load(fullfile(folder,[filename,'_regobj.mat']))
    fig = fopen(fullfile(folder,imgfolderlist(i).name),'r');
    [Ly, Lx] = size(RegPara.mimg);
    clear mov; clear Region_full; clear trace_full;
    mov = fread(fig,Ly*Lx*metadata.numframe,'*uint8');
    mov = reshape(mov, Ly,Lx,metadata.numframe);
    fclose(fig);
    
    % gaussian filtering
    tic
    movF = MovGaussFilter(mov, G, subsample);
    toc
    %%%% find neuron ROIs using kurtosis methods
    tic
    [movKurt_ori, rois, regions] = batchKurt_analysis(movF, minArea, displayit, SegmentPara);
    toc
    clear movF; 
    %%%% remove ROI partially out of imaging window during time lapse based on registration results
    if ~isfield(RegPara,'yrange')
        RegPara.yrange = 1:Lx;
        RegPara.xrange = 1:Ly;
    end
    
    ROIrange = zeros(Ly,Lx);
    ROIrange(:,RegPara.yrange)=ROIrange(:,RegPara.yrange)+1;
    ROIrange(RegPara.xrange,:)=ROIrange(RegPara.xrange,:)+1;
    ROIrange(ROIrange<2)=0;
    regions1 = regions;
    regions1(ROIrange==2) = 0;
    id = unique(regions1(:));
    id(id==0) = [];
    if ~isempty(id)
        rois(:,:,id) = [];
        for j = 1:length(id)
            regions(regions==id(j)) = 0;
        end
    end
    Region_full.movKurt = movKurt_ori;
    Region_full.rois = rois;
    Region_full.regions = regions;
    Region_full.G = G;
    Region_full.minArea = minArea;
    Region_full.SegmentPara = SegmentPara;
    save([filename,'_traces.mat'],'Region_full');
    
    rois = Region_full.rois;
    tic
    [tcraw, tc, tcraw_pca, tc_pcasub] = kyoGetTCs_sutract_pca(mov, rois, subtractPara);
    toc   
    Traces_full.tcraw = tcraw;
    Traces_full.tc = tc;
    Traces_full.tcraw_pca = tcraw_pca;
    Traces_full.tc_pcasub = tc_pcasub;
    Traces_full.subtractPara = subtractPara;  
    save([filename,'_traces.mat'],'Traces_full','-append');
end 

%% generate ROI using suite2p
subtractPara = [2,5]; % neuropil correction parameters
imgfolderlist = dir(fullfile(imgfolder,'*_reg.bin'));
%
SVDROIPara.useGPU                 = 0; % if you can use a GPU in matlab this accelerate registration approx 3 times
SVDROIPara.sig                    = 0.8;  % spatial smoothing length for clustering; encourages localized clusters
SVDROIPara.nSVDforROI             = 500; % Number of principal components to keep
SVDROIPara.NavgFramesSVD          = 4000; % how many (binned) timepoints to do the SVD based on
SVDROIPara.niterclustering        = 30;   % how many iterations of clustering
SVDROIPara.ShowCellMap            = 1;
SVDROIPara.Nk0                    = 800;  % how many clusters to start with
SVDROIPara.Nk                     = 100;  % how many clusters to end with
SVDROIPara.resultsavepath = datafolder;
SVDROIPara.minarea = 100; % minimal cell size (was 50)
SVDROIPara.maxarea = 3000; % maximal cell size (was 500)
% RegPara.useGPU                 = 0; % if you can use a GPU in matlab this accelerate registration approx 3 times
% RegPara.sig                    = 0.8;  % spatial smoothing length for clustering; encourages localized clusters
% RegPara.nSVDforROI             = 1000; % Number of principal components to keep
% RegPara.NavgFramesSVD          = 4000; % how many (binned) timepoints to do the SVD based on
% RegPara.niterclustering        = 50;   % how many iterations of clustering
% RegPara.ShowCellMap            = 1;
% RegPara.Nk0                    = 800;  % how many clusters to start with
% RegPara.Nk                     = 400;  % how many clusters to end with
% RegPara.resultsavepath = datafolder;
% RegPara.minarea = 100; % minimal cell size (was 50)
% RegPara.maxarea = 1600; % maximal cell size (was 500)
for i = 1:length(imgfolderlist)
    i
    filename = strtrim(imgfolderlist(i).name(1:end-8))
    load(fullfile(imgfolder,[filename,'_regobj.mat']))
    clear U; clear Sv; clear V
    tic
    [SVDROIPara, U, Sv] = get_svdForROI_v2(SVDROIPara, RegPara, metadata);
%     [RegPara, U, Sv] = get_svdForROI_v2(RegPara, 'CNMF');
    [SVDROIPara, stat, res] = fast_cluster_neuropil_v2(SVDROIPara, RegPara, U, Sv);
%     [RegPara, stat, res] = fast_clustering_with_neuropil2(RegPara, U, Sv);
    toc
    %%%% ROI based on morphology selection
    roi = []; k=1;
    for ii = 1:length(stat)
        roi1 = zeros(size(RegPara.mimg));
        roi1(stat(ii).ipix)=1;
        bw = bwmorph(roi1,'clean');
        bw = bwlabel(bw);
        if max(bw(:))>0
            for j = 1:max(bw(:))
                [a,b] = find(bw==j);
                r1 = length(find(bw==j));
                if r1>SVDROIPara.minarea && r1<SVDROIPara.maxarea
                    roi2 = zeros(size(Reg.mimg));
                    roi2(bw==j)=1;
                    roi(:,:,k) = logical(roi2);
                    k = k+1;            
                end
            end
        end
    end
    roi_suite = logical(roi);
    save([filename,'_SVDtraces.mat'],'U','Sv','res','stat','SVDROIPara','roi_suite');
    %%%% extract calcium traces
    tic
    [tcraw, tc, tcraw_pca, tc_pcasub] = kyoGetTCs_sutract_pca(mov, roi_suite, subtractPara);
    toc 
    clear mov
    Traces_full.tcraw = tcraw;
    Traces_full.tc = tc;
    Traces_full.tcraw_pca = tcraw_pca;
    Traces_full.tc_pcasub = tc_pcasub;
    Traces_full.subtractPara = subtractPara;  
    save([filename,'_traces.mat'],'roi_suite','Traces_full','-append');
end