% Example Script for PSF Deconvolution

% Load EPI testdata
load('sampleData.mat');
kspData = data.kspData;
fieldMap = data.fieldMap;

% Set parameters
params.nx = size(kspData,1); params.ny = size(kspData,2); 
params.Te = 0.034; params.Ta = 0.15; 
params.yres = 0.192*1e3/params.nx;
params.alpha = 20;
scan.type = 'ge-epi';
t2star = 20;

% Reconstruct data
imageDistorted     = zeros(size(kspData));
reconImageZeroFill = zeros(size(kspData));
imageConjFill      = zeros(size(kspData));
reconImageConjFill = zeros(size(kspData));

w = waitbar(0,'Correcting images ...');
for Nslice = 1:size(kspData,3)
    waitbar(Nslice/size(kspData,3), w, 'Correcting images ...');
    
    % (A) Zerofilling
    % Transform k-space data to image domain
    imageDistorted(:,:,Nslice) = fftshift(fftn(ifftshift(kspData(:,:,Nslice))));
    scan.pf = 0.56; scan.pftype = 'zerofill'; scan.direction = 'down';
    reconImageZeroFill(:,:,Nslice) = PSF_correction_method_patzig(imageDistorted(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
    
    % (B) Conjugatefilling
    % ksp_tmp = ifftshift(ifftn(fftshift(imageEPI(:,:,Nslice))));
    ksp_tmp = kspData(:,:,Nslice);
    kspConj = rot90(conj(kspData(:,:,Nslice)),2);
    if strcmp(scan.direction, 'up')
        ksp_tmp(ceil(scan.pf*size(ksp_tmp,1)):end,:) = kspConj(ceil(scan.pf*size(ksp_tmp,1)):end,:);
    elseif strcmp(scan.direction, 'down')
        ksp_tmp(1:end-ceil(scan.pf*size(ksp_tmp,1)),:) = kspConj(1:end-ceil(scan.pf*size(ksp_tmp,1)),:);
    end
    imageConjFill(:,:,Nslice) = fftshift(fftn(ifftshift(ksp_tmp)));
    
    scan.pftype = 'conjfill';
    reconImageConjFill(:,:,Nslice) = PSF_correction_method_patzig(imageConjFill(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
end
delete(w);

% Plot Results
figure
imagesc([abs(imageDistorted(:,:,Nslice))/max(max(abs(imageDistorted(:,:,Nslice)))) ...
         abs(reconImageZeroFill(:,:,Nslice))/max(max(abs(reconImageZeroFill(:,:,Nslice)))) ...
         abs(reconImageConjFill(:,:,Nslice))/max(max(abs(reconImageConjFill(:,:,Nslice))))])
title(['Slice ' num2str(Nslice) ': Distorted image - Corrected image (zero filling) - Corrected image (conjugate filling)'])
colormap(gray); axis image; axis off;

%% Example for DEPICTING (not included in the testdata)
params.nx = size(imageDEPICTING,1); params.ny = size(imageDEPICTING,2);
params.Te = 0.005; params.Ta = 0.15; 
params.yres = 0.192*1e3/params.nx;
scan.type = 'depicting';
t2star = 0.2;

Nslice = 10;
reconImage{Nslice} = PSF_correction_method_patzig(imageDEPICTING(:,:,Nslice), fieldMap(:,:,Nslice), t2star*ones(size(fieldMap(:,:,Nslice))), params, scan);
