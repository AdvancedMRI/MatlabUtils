function [x_1,x_2,x] = nii4ensight(stem,scale,shift)
% convert .mat file into a nifti file for use with ensight
% import nii file to ensight, e.g. edit ensight_template or other saved ensight context to point
% to nii file
% alex beckett, 2020

x_1 = [];
x_2 = [];
% code to convert raw 4dflow data into an nii file suitable for analysis in
% enight
if nargin <1
    disp('Specify a MID folder with a measyaps_protocol.evp file')
    % other option - open gui for dir selection
    % create default state if no evp file exists?
    return
end

if nargin<2 || isempty(scale)
    scale=1;
    % scaling factor as needed to work with prev analyses
    % might be less important, consider moving down
end

% make record of slice orientation for checking we are using right code 
cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.asSlice\[0].sNormal',char(39),' | awk ',char(39),'{print $1}',char(39)];
[~,sDir] = system(cmd);
sDir = sDir(end-3:end-1);

% find out slice orientation from .evp
if strcmp(sDir,'Tra')
    disp('Looks like this this is Axial data, continuing as planned')
    % add some stuff to preplan dim swaps and phase resigning here -
    % hardcode for now
    dimSwap = [2 1 3 4 5];
    phaseFlip = [1 1 1];
elseif strcmp(sDir,'Sag')
    disp('Looks like this this is Sagittal data, continuing as planned')
    % add some stuff to preplan dim swaps and phase resigning here
    dimSwap = [3 2 1 4 5];
    phaseFlip = [1 1 -1];
else
   disp('I will write code that copes with whatever this is one day! proceed with caution!')
   phaseFlip = [1 1 1];
end

% load data from .mat file
x=load([stem,'/result3']);

% get relevant data from .evp file
cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.asSlice\[0].dThickness',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,dThick]=system(cmd);
    
cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.asSlice\[0].dPhaseFOV',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,dPhaseFOV]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'alTR\[0]',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,TR]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'alTE\[0]',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,TE]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sKSpace.lPhaseEncodingLines  ',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,lPhaseEncodingLines]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sKSpace.lPartitions  ',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,lPartitions]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sAngio.sFlowArray.asElm\[0].nVelocity  ',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,venc]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.lSize  ',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,slabs]=system(cmd);

cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sGroupArray.asGroup\[0].dDistFact   ',char(39),' | awk ',char(39),'{print $3}',char(39)];
[~,gap]=system(cmd);
gap = str2double(gap);

% slab order comparing MB1 and MB2 is ambiguos - save this line from the
% header into the description of the nii file to aid interpretation
cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.asSlice\[0].sPosition.d',sDir,char(39)];
[~,note] = system(cmd);
note=note(1:end-1);

% gap and slice shift will interact in a weird way, look into it
if nargin<3  
    shift = gap;
end

% calculate resolution and fov from evp
inPlane=str2double(dPhaseFOV)/str2double(lPhaseEncodingLines);
throughPlane = str2double(dThick)/str2double(lPartitions);
slabs = str2double(slabs)-1;
slices = str2double(lPartitions);

voxel=[inPlane inPlane throughPlane]; % mm 
venc=str2double(venc)*10; % convert to mm/s;
tr=str2double(TR)/10^6; %convert to s

% rescale mag if requested
x.mag_all=x.mag_all*scale;

% convert phase to velocity - need a way to evaluate which phase dirs are which?
x.phase_all=x.phase_all/pi*venc;
for i=1:3
x.phase_all(:,:,:,:,i)=x.phase_all(:,:,:,:,i)*phaseFlip(i);
end

% swap 3rd and 4th dim so time is the 4th dimension
x.mag_all=(permute(x.mag_all,[1 2 4 3]));
x.phase_all=(permute(x.phase_all,[1 2 4 3 5]));

% slab 1 will always be the first n partitions
x_1.mag_all=x.mag_all(:,:,1:slices,:);
x_1.phase_all=x.phase_all(:,:,1:slices,:,:);

% calculate absolute velocity r for masking
[x_1.th,x_1.phi,x_1.r]=cart2sph(x_1.phase_all(:,:,:,:,2),x_1.phase_all(:,:,:,:,3),x_1.phase_all(:,:,:,:,1));

% create parameter maps for masking
magMean=mean(x_1.mag_all,4);
veloStd=std(x_1.r,[],4);

% currently noise values for masking are hardcoded based on Bock et al -
% consider parametizing

% noise masking via 9 +/- 3 max mag (Bock)
magMask=magMean>(max(magMean)*.09);

% noise masking via 20 +/- 4 max velo std (Bock)
% do sep for each component (a la Walker)?
veloMask=veloStd<(max(veloStd(:))*.3);

% combine noise masks
allMask=magMask.*veloMask;

% stationary tissue masks
% via 10 +/-2 max velo std (Bock)
statMask=veloStd<(max(veloStd(:))*.05);

% via 15th percentils of velo std (walker)
% should be per component, try later
statMask_walker=veloStd<prctile(veloStd(veloStd>0),15);

% fit plane to each component
for i=1:3
    % this is wrong, need to either do each frame, or pick a single frame?
    mP(:,:,:,i)=mean(x_1.phase_all(:,:,:,:,i),4);
    thisSP=std(x_1.phase_all(:,:,:,:,i),[],4);
    sP(:,:,:,i)=thisSP;
    % stationary mask for each phase?
    sPm(:,:,:,i)=sP(:,:,:,i)<(max(thisSP(:))*.1);
    mPv=mP(sPm(:,:,:,i));
    %[r,c,v]=ind2sub(size(statMask),find(statMask));
    [r,c,v]=ind2sub(size(sPm(:,:,:,i)),find(sPm(:,:,:,i)));
    DM = [r c v ones(length(r),1)];
    B=DM\mPv;
    [X,Y,Z]=ndgrid(1:size(x_1.phase_all,1),1:size(x_1.phase_all,2),1:size(x_1.phase_all,3));
    x_1.P(:,:,:,1,i) = B(1)*X + B(2)*Y + B(3)*Z +B(4)*ones(size(X));
    x_1.mP=mP;
end
% subtract out phase bge calculated above
x_1.phase_all_corr=x_1.phase_all-repmat(x_1.P,[1 1 1 size(x_1.phase_all,4) 1]);

% mask out noise
x_1.phase_all=x_1.phase_all.*repmat(allMask,[1 1 1 size(x_1.phase_all,4) 3]);
x_1.mag_all=x_1.mag_all.*repmat(allMask,[1 1 1 size(x_1.phase_all,4)]);

x_1.phase_all_corr=x_1.phase_all_corr.*repmat(allMask,[1 1 1 size(x_1.phase_all,4) 3]);

% mask out stationary (doesnt work great)
x_1.phase_all_stat=x_1.phase_all_corr.*repmat(~statMask,[1 1 1 size(x_1.phase_all,4) 3]);

% Ipc MRA
% recalculate absolute velocity for MRA
[x_1.th,x_1.phi,x_1.r]=cart2sph(x_1.phase_all(:,:,:,:,2),x_1.phase_all(:,:,:,:,3),x_1.phase_all(:,:,:,:,1));
[x_1.th_corr,x_1.phi_corr,x_1.r_corr]=cart2sph(x_1.phase_all_corr(:,:,:,:,2),x_1.phase_all_corr(:,:,:,:,3),x_1.phase_all_corr(:,:,:,:,1));

% multiply by image intensity (IpcMRA) and self weight (sum in ensight)
x_1.IpcMRAw=(x_1.mag_all.*x_1.r).^2;
x_1.IpcMRAw_corr=(x_1.mag_all.*x_1.r_corr).^2;

% display results
figure('Name',stem)
for i=1:3
    subplot(3,4,i),imagesc((squeeze(x_1.phase_all(:,:,end/2,end,i).*allMask(:,:,end/2))),[-venc/2 venc/2]),axis('image'),colormap('gray')
    subplot(3,4,i+4),imagesc((squeeze(x_1.P(:,:,end/2,i))),[-venc/2 venc/2]),axis('image'),colormap('gray')
    subplot(3,4,i+8),imagesc((squeeze(x_1.phase_all_corr(:,:,end/2,end,i).*allMask(:,:,end/2))),[-venc/2 venc/2]),axis('image'),colormap('gray')
end

subplot(3,4,12),imagesc((squeeze(max(mean(x_1.IpcMRAw_corr,4),[],3))),[0 .05]),axis('image'),colormap('gray')


% now process second slab if this is MB2 data

if slabs %split MB into 2 - sometimes requires an offset in the second slab (make an argument?)
    disp('Processing second slab')
    if shift <1
        x_2.mag_all=x.mag_all(:,:,[slices+(slices*shift)+1:slices*2,slices+1:slices+(slices*shift)],:); % this needs un hard  coding - calculate based on gap? add a parameter?
        x_2.phase_all=x.phase_all(:,:,[slices+(slices*shift)+1:slices*2,slices+1:slices+(slices*shift)],:,:);
    else
        x_2.mag_all=x.mag_all(:,:,[slices+1:slices*2],:); % this needs un hard  coding - calculate based on gap? add a parameter?
        x_2.phase_all=x.phase_all(:,:,[slices+1:slices*2],:,:);
    end
    [x_2.th,x_2.phi,x_2.r]=cart2sph(x_2.phase_all(:,:,:,:,2),x_2.phase_all(:,:,:,:,3),x_2.phase_all(:,:,:,:,1));
    
    magMean=mean(x_2.mag_all,4);
    magMask=(magMean>prctile(magMean(magMean>0),9));
    veloStd=std(x_2.r,[],4);
    veloMask=veloStd<(max(veloStd(:))*.3);
    allMask=magMask.*veloMask;
    x_2.phase_all=x_2.phase_all.*repmat(allMask,[1 1 1 size(x_2.phase_all,4) 3]);
    x_2.mag_all=x_2.mag_all.*repmat(allMask,[1 1 1 size(x_2.phase_all,4)]);
    
    [x_2.th,x_2.phi,x_2.r]=cart2sph(x_2.phase_all(:,:,:,:,2),x_2.phase_all(:,:,:,:,3),x_2.phase_all(:,:,:,:,1));
    x_2.IpcMRAw=(x_2.mag_all.*x_2.r).^2;
    
    if 1 %bg phase corr -always do, argument? not sure
        for i=1:3
            % this is wrong, need to either do each frame, or pick a single frame?
            mP(:,:,:,i)=mean(x_2.phase_all(:,:,:,:,i),4);
            thisSP=std(x_2.phase_all(:,:,:,:,i),[],4);
            sP(:,:,:,i)=thisSP;
            % stationary mask for each phase?
            sPm(:,:,:,i)=sP(:,:,:,i)<(max(thisSP(:))*.1);
            mPv=mP(sPm(:,:,:,i));
            %[r,c,v]=ind2sub(size(statMask),find(statMask));
            [r,c,v]=ind2sub(size(sPm(:,:,:,i)),find(sPm(:,:,:,i)));
            DM = [r c v ones(length(r),1)];
            B=DM\mPv;
            [X,Y,Z]=ndgrid(1:size(x_2.phase_all,1),1:size(x_2.phase_all,2),1:size(x_2.phase_all,3));
            x_2.P(:,:,:,1,i) = B(1)*X + B(2)*Y + B(3)*Z +B(4)*ones(size(X));
            x_2.mP=mP;
        end
        x_2.phase_all_corr=x_2.phase_all-repmat(x_2.P,[1 1 1 size(x_2.phase_all,4) 1]);
        x_2.phase_all_corr=x_2.phase_all_corr.*repmat(magMask,[1 1 1 size(x_2.phase_all,4) 3]);
        [x_2.th_corr,x_2.phi_corr,x_2.r_corr]=cart2sph(x_2.phase_all(:,:,:,:,2),x_2.phase_all(:,:,:,:,3),x_2.phase_all(:,:,:,:,1));
        x_2.IpcMRAw_corr=(x_2.mag_all.*x_2.r_corr).^2;

        
        
    end
    
    % display results
    figure('Name',stem)
    for i=1:3
        subplot(3,4,i),imagesc((squeeze(x_2.phase_all(:,:,end/2,end,i).*allMask(:,:,end/2))),[-venc/2 venc/2]),axis('image'),colormap('gray')
        subplot(3,4,i+4),imagesc((squeeze(x_2.P(:,:,end/2,i))),[-venc/2 venc/2]),axis('image'),colormap('gray')
        subplot(3,4,i+8),imagesc((squeeze(x_2.phase_all_corr(:,:,end/2,end,i).*allMask(:,:,end/2))),[-venc/2 venc/2]),axis('image'),colormap('gray')
    end
    
    % recombine MB
    gapSize = size(x_1.phase_all);
    gapSize(3) = gapSize(3)*ceil(gap);
    x.phase_all=cat(3,x_1.phase_all,zeros(gapSize),x_2.phase_all); % base zeros on gap factor
    x.mag_all=cat(3,x_1.mag_all,zeros(gapSize(1:4)),x_2.mag_all);
    x.IpcMRAw=cat(3,x_1.IpcMRAw,zeros(gapSize(1:4)),x_2.IpcMRAw);
    
    % swap dims depending on slice orientation
    
    % slab order comparing MB1 and MB2 is ambiguos - save this line from the
    % header into the description of the nii file to aid interpretation
    cmd=['cat ',stem,'/measyaps_protocol.evp | grep -a ',char(39),'sSliceArray.asSlice\[1].sPosition.d',sDir,char(39)];
    [~,note2] = system(cmd);
    note2=note2(1:end-1);
    % save
    
    vol = cat(5,x.mag_all,x.phase_all,x.IpcMRAw);
    switch sDir
        case 'Tra'
            vol=permute(vol,[2 1 3 4 5]);
            voxel = voxel([2 1 3]);
        case 'Sag'
            vol=flip(permute(vol,[3 2 1 4 5]),3);
            voxel = voxel([3 2 1]);
        otherwise
    end
    
    nii=make_nii(vol,voxel);
    nii.hdr.dime.pixdim(5)=tr;
    nii.hdr.dime.xyzt_units = 10;
    nii.hdr.hist.descrip=note;
    save_nii(nii,[stem,'_',sDir,'_MB',num2str(slabs+1),'_slab1+2.nii']);
    
    vol = cat(5,x_2.mag_all,x_2.phase_all_corr,x_2.IpcMRAw_corr);
    switch sDir
        case 'Tra'
            vol=permute(vol,[2 1 3 4 5]);
            voxel = voxel([2 1 3]);
        case 'Sag'
            vol=flip(permute(vol,[3 2 1 4 5]),3);
            voxel = voxel([3 2 1]);
        otherwise
    end
    
    nii=make_nii(vol,voxel);
        nii.hdr.dime.pixdim(5)=tr;
        nii.hdr.dime.xyzt_units = 10;
        save_nii(nii,[stem,'_',sDir,'_MB',num2str(slabs+1),'_slab2_corr.nii']);
    
    vol = cat(5,x_2.mag_all,x_2.phase_all,x_2.IpcMRAw);
    switch sDir
        case 'Tra'
            vol=permute(vol,[2 1 3 4 5]);
            voxel = voxel([2 1 3]);
        case 'Sag'
            vol=flip(permute(vol,[3 2 1 4 5]),3);
            voxel = voxel([3 2 1]);
        otherwise
    end
    
	nii=make_nii(vol,voxel);
    nii.hdr.dime.pixdim(5)=tr;
    nii.hdr.dime.xyzt_units = 10;
    nii.hdr.hist.descrip=note2;
    save_nii(nii,[stem,'_',sDir,'_MB',num2str(slabs+1),'_slab2.nii']);
end

% finalize slab 1
% swap dims depending on slice orientation

vol = cat(5,x_1.mag_all,x_1.phase_all,x_1.IpcMRAw);
switch sDir
    case 'Tra'
        vol=permute(vol,[2 1 3 4 5]);
        voxel = voxel([2 1 3]);
    case 'Sag'
        vol=flip(permute(vol,[3 2 1 4 5]),3);
        voxel = voxel([3 2 1]);
    otherwise
end

nii=make_nii(vol,voxel);
nii.hdr.dime.pixdim(5)=tr;
nii.hdr.dime.xyzt_units = 10;
nii.hdr.hist.descrip=note;
save_nii(nii,[stem,'_',sDir,'_MB',num2str(slabs+1),'_slab1.nii']);


vol = cat(5,x_1.mag_all,x_1.phase_all_corr,x_1.IpcMRAw_corr);
switch sDir
    case 'Tra'
        disp('Swapping dims for Axial slices')
        vol=permute(vol,[2 1 3 4 5]);
        voxel = voxel([2 1 3]);
    case 'Sag'
        disp('Swapping dims for Sagittal slices')
        vol=flip(permute(vol,[3 2 1 4 5]),3);
        voxel = voxel([3 2 1]);
    otherwise
        disp('Unknown slices')
end
nii=make_nii(vol,voxel);
nii.hdr.dime.pixdim(5)=tr;
nii.hdr.dime.xyzt_units = 10;
nii.hdr.hist.descrip=note;
save_nii(nii,[stem,'_',sDir,'_MB',num2str(slabs+1),'_slab1_corr.nii']);
