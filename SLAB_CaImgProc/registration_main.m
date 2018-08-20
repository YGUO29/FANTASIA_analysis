% registration modified from suite2p
folder = 'C:\Users\slab\Desktop\CaImgProc\testdata';
imgfolderlist = dir(fullfile(folder,'*.tif'));
metadata.expdate = ' '; % experiment date
metadata.stim = 'MovingSquare'; % visual stimulus
metadata.Regsavepath = 'C:\Users\slab\Desktop\CaImgProc\testdata\processed'; % path for saving registration results
RegPara.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times
RegPara.NimgFirstRegistration  = 100;
RegPara.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
RegPara.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
RegPara.NiterPrealign          = 20;
RegPara.RawPrecision           = 'uint8';

for i = 1:size(imgfolderlist,1)
        filename = imgfolderlist(i).name
        imglist = imfinfo(fullfile(folder,filename));
        data1 = imread(fullfile(folder,filename),1);
        metadata.imgsize = size(data1);
        metadata.numframe = length(imglist);
        metadata.name = strtrim(filename);
        % initial alignment 
        clear data
        f1 = round(rand(1)*length(imglist)/2);
        for j = f1+1:f1+RegPara.NimgFirstRegistration
            data(:,:,j-f1) = imread(fullfile(folder,filename),j);
        end
        data = single(data);
        RegPara = align_iterative_modified
        (data, RegPara);           
        clear data
        
        % register whole movie
        RegPara.RegFile = fullfile(metadata.Regsavepath, [metadata.name(1:end-14),'_reg.bin']);
        regdir = fileparts(RegPara.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        tic
        ix0 = 0; Nbatch = 500; 
        dreg = zeros([size(RegPara.mimg), length(imglist)],class(data1));
        
        dsall = []; CorrAll = [];
        while ix0<length(imglist)
            clear mov            
            indxr = ix0 + (1:Nbatch);            
            indxr(indxr>length(imglist)) = [];
            mov = zeros([size(data1), length(indxr)],'uint8');
            for j = indxr(1):indxr(end)
                mov(:,:,j-ix0) =  imread(fullfile(folder,filename),j);
            end
            
            [ds, Corr]  = registration_offsets_modified(mov, RegPara, 0);
            CorrAll = cat(1,CorrAll, Corr);
            dsall = cat(1,dsall, ds);

            dreg(:, :, indxr) = register_movie(mov, RegPara, ds);
            ix0 = ix0 + Nbatch;
            indxr(end)
        end
        toc
        
        RegPara.dsall = dsall;
        RegPara.CorrAll = CorrAll;
        RegPara = effectiveRange(RegPara);
        % save result
        clear fid
        fid = fopen(RegPara.RegFile, 'w');
        fwrite(fid, dreg, class(data1));
        fclose(fid);            
        % save reg sample
        sampledreg = dreg(:,:,1:15:end);
        clear mov
        imwrite(sampledreg(:,:,1), fullfile(metadata.Regsavepath, [metadata.name,'_RegSample.tif']))
        for k = 2:size(sampledreg,3)
            imwrite(sampledreg(:,:,k), fullfile(metadata.Regsavepath, [metadata.name,'_RegSample.tif']), 'writemode', 'append');
        end
        clear dreg; clear sampledreg
        save(fullfile(metadata.Regsavepath, [metadata.name(1:end-14),'_regobj.mat']),'metadata','RegPara');
toc
end

%% read bin data
fid = fopen(RegPara.RegFile, 'r');
data = fread(fid, metadata.imgsize(2)*metadata.imgsize(1)*metadata.numframe, '*uint8');
mov = reshape(data, metadata.imgsize(1),metadata.imgsize(2),metadata.numframe);