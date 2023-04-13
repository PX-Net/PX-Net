function [ V VS ] = dicomreadfolder2(folder)
%DICOMREADVOLUME   Read a zip-file with DICOM slices as volume.
%   [V VS] = DICOMREADVOLUME(FNAME) loads the 3-D volume data V from
%   the volume dicom slices in the zip-file FNAME. 
%   The voxel size is returned by VS.
%
%   See also DICOMWRITEVOLUME, DICOMWRITE, DICOMREAD
%
%   Author: medical-image-processing.com


% extract the zip file to the temporary directory

dcm_list=[folder,'/*.dcm'];
dcm_names=dir(dcm_list);
N = numel(dcm_names);
if N==1   
    dcm_names(1).name=[folder,'/',dcm_names(1).name];
    I1 = dicominfo(dcm_names(1).name);
    Vt=squeeze(dicomread(dcm_names(1).name));
    if length(size(Vt))>=3
        V=Vt;
    end
    ps=I1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    st=I1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
    VS=[ps' st];
end

for i=1:N
    dcm_names(i).name=[folder,'/',dcm_names(i).name];
end


% read first slice for determining slice properties
S = squeeze(dicomread(dcm_names(1).name));
I1 = dicominfo(dcm_names(1).name);

% voxel size information
if isempty(I1.PixelSpacing(1))
    I1.PixelSpacing(:)=[1;1];
end
if isempty(I1.SliceThickness)
    I1.SliceThickness=1;
end
if isempty(I1.RescaleIntercept(:))
    I1.RescaleIntercept(:)=0;
end

VS = [I1.PixelSpacing(:) ; I1.SliceThickness ; I1.RescaleIntercept(:)];

% slice size and datatype
sz = size(S);
tp = class(S);

% pre-allocate data
VT = zeros([sz N], tp);
V = VT;
POS = zeros(N,2);

% load each slice and its properties
for i=1:N
    VT(:,:,i) = squeeze(dicomread(dcm_names(i).name));
    info = dicominfo(dcm_names(i).name);
    if isfield(info, 'ImagePositionPatient')
        POS(i,:) = [info.ImagePositionPatient(3) i];
    else
        POS(i,:) = [i i];
    end
%     delete(fnames{i});
end

% resort the slices according to the image position
POS = sortrows(POS,-1);
for i=1:N
   V(:,:,i) = VT(:,:,POS(i,2)); 
end


