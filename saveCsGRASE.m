function saveCsGRASE(stem,voxSize,qflag,zflag,pflag)
% function to convert recon'd GRASE data in .mat format into .nii for
% further processing
% requires MATLAB Nifti toolbox (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
% requires mrTools (https://gru.stanford.edu/doku.php/mrtools/download)
% stem = stem of the files you with to convert (e.g. 'MID*')
% voxSize = 4 element array of voxel size in [mm mm mm sec]
% qflag = flag to create a qform, required for mrTools but causes issues for FSL
% zflag = flip image in partition phase encode direction during creation
% plag = flip image in inplane phase encode direction during creation

if nargin<1 || isempty(stem)
    stem='*.mat';
end

% set default voxel voxSize to 1mm iso
if nargin<2
    voxSize=[1 1 1 1];
end

% by default do not set-up an MLR valid qform
if nargin<3
    qflag=0;
end

% by default do not flip z axis
if nargin<4
    zflag=0;
end

% by default we flip pe axis
if nargin<4
    pflag=1;
end

d=dir(stem);
for i=1:length(d)
    load(d(i).name)
    s=whos('img*');
    tmp=eval(s.name);
    tmp = (permute(tmp(:,:,:,:),[2 1 3 4])); %used to have to flip in slice di, now we flip in short PE dir
    if pflag
        %flip slice axis if needed
        tmp = tmp(:,end:-1:1,:,:);
    end
    if zflag
        %flip slice axis if needed
        tmp = tmp(:,:,end:-1:1,:);
    end
    tmp = (tmp/max(tmp(:)))*4095;
    nii=make_nii(tmp,voxSize(1:3),[],4);
    nii.hdr.dime.pixdim(5)=voxSize(4);
    save_nii(nii,[stripext(d(i).name),'.nii']);
    clear img*
end

% optional step to set-up a valid qform, needed for MLR tools
% may cause issues for FSL type analyses

if qflag
    d=dir('*.nii');
    for i=1:length(d)
    [a,b]=cbiReadNifti(d(i).name);
    b.xyzt_units=10;
    b.qform_code=1;
    cbiWriteNifti(b.hdr_name,a,b);
    end
end
