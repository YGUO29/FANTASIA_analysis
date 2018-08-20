function [ops, U, Sv, V] = get_svdForROI_v2(ops, RegPara, matadata)
Ly = size(RegPara.mimg,1);
Lx = size(RegPara.mimg,2);
ntotframes = matadata.numframe;
ops.NavgFramesSVD  = min(ops.NavgFramesSVD, ntotframes);
nt0 = ceil(ntotframes / ops.NavgFramesSVD); % number of frames to average
ops.NavgFramesSVD = floor(ntotframes/nt0);

fid = fopen(RegPara.RegFile, 'r');
data = fread(fid,  Ly*Lx*ntotframes, '*uint8');
data = reshape(data, Ly, Lx, []);
irange = 1:nt0*floor(size(data,3)/nt0);
if irange(end)<ntotframes
    data(:,:,irange(end)+1:ntotframes) = [];
end

tic
mov = zeros(Ly, Lx, ops.NavgFramesSVD, 'single');
ix = 0;
while ix<irange(end)/nt0
    if (ix+100)*nt0<=irange(end) % batch size 100 (compute 100 average image at a time)
        data1 = reshape(data(:,:,ix*nt0+1:(ix+100)*nt0), Ly, Lx, nt0, []);
    else
        data1 = reshape(data(:,:,ix*nt0+1:end), Ly, Lx, nt0, []);
    end
        davg = single(squeeze(mean(data1,3)));
        clear data1
        davg = davg - repmat(mean(davg,3), 1, 1, size(davg,3));
        mov(:,:,ix + (1:size(davg,3))) = davg;
        ix = ix + size(davg,3);
        clear davg            
end
fclose(fid);
clear data
toc
%%
mov(:, :, (ix+1):end) = [];
% if exist(RegPara.yrange)
%     mov = mov(RegPara.yrange, RegPara.xrange, :);
% end
[c,r,~] = size(mov);
%% SVD options
ops.nSVDforROI = min(ops.nSVDforROI, size(mov,3));

if ops.sig>0.05
	for i = 1:size(mov,3)
	   I = mov(:,:,i);
	   I = my_conv(my_conv(I',ops.sig)', ops.sig);
	   mov(:,:,i) = single(I);
	end
end

mov = reshape(mov, [], size(mov,3));
mov = mov./repmat(mean(mov.^2,2).^.5, 1, size(mov,2));
COV = mov' * mov/size(mov,1);
ops.nSVDforROI = min(size(COV,1)-2, ops.nSVDforROI);

if ops.useGPU && size(COV,1)<1.2e4
    reset(gpuDevice); 
    g = gpuDevice; 
    disp(g.FreeMemory);
    [V, Sv, ~]      = svd(gpuArray(double(COV)));
    V               = single(V(:, 1:ops.nSVDforROI));
    Sv              = single(diag(Sv));
    Sv              = Sv(1:ops.nSVDforROI);
    %
     Sv = gather(Sv);
     V = gather(V);
    reset(gpuDevice); 
    g = gpuDevice; 
%     disp(g.FreeMemory);
else
    [V, Sv]          = eigs(double(COV), ops.nSVDforROI);
    Sv              = single(diag(Sv));
end

U               = normc(mov * V);
U               = single(U);
U = reshape(U,c,r,[]);

