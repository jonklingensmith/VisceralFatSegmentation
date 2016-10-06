%------------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral researcher
% Institution: CNG at TTU in a working association with TTNI
% Date: 1 JUL 2015
% Updated: 31 JAN 2016
%------------------------------------------------------------------------

clear all, close all;

addpath(genpath(pwd))

% subject: MF0201
% session: 1012151
% date:    12 OCT 2015

% image generation control
showProgress = false; %true;
batchShow    = ~showProgress;

% NOTE: For batch runs over many slices should set showProgress to 'false'
%       For examining the analysis of an individual slice,
%       set showProgress to 'true' to see more infomation.
%       For individual slice run comment out the for i=1:Nslices line and
%       its 'end' line. Also comment out the saving of the .mat file.
%       Slices at the bottom of the heart cannot be viewed individually,
%       since they are masked from the slice above. Must run the for loop
%       from heartMax to the slice of interest.

%% scan info

NIFTI_file_name_F_upper = '20151012_143006t1vibedixontrap4bh320s004a1001.nii.gz';
NIFTI_file_name_W_upper = '20151012_143006t1vibedixontrap4bh320s005a1001.nii.gz';
NIFTI_file_name_F_lower = '20151012_143006t1vibedixontrap4bh320s009a1001.nii.gz';
NIFTI_file_name_W_lower = '20151012_143006t1vibedixontrap4bh320s010a1001.nii.gz';

voxelSize_upper = [1.40625, 1.40625, 2.5];
interslice_spacing_fraction = 0.2;
voxelVolume_upper = prod(voxelSize_upper)*(1+interslice_spacing_fraction);
voxelSize_lower = [1.40625, 1.40625, 2.5];
voxelVolume_lower = prod(voxelSize_upper)*(1+interslice_spacing_fraction);
X_shift = -round(2.2/voxelSize_lower(1));
Y_shift = -round(9.1/voxelSize_lower(2));

Z_L1_L2u    = 11; % middle of L1-L2 intervertebral disk slice (upper scan)
%Z_L1        = 13;  % bottom of L1 vertebra slice
Z_T4_T5t       = 112; % top of T4-T5 intervertebral disk slice
upperSlices = Z_T4_T5t:-1:Z_L1_L2u;
Z_T8_T9d       = 55;  % top of T9 vertebra disk slice (use for bottom of thoracic cavity) 
% NOTE: diaphraim is at slices = 
% - register from the bottom of T8
heartMax    = Z_T4_T5t-70+1;
heartBottom = Z_T4_T5t-59+1;
diaphram    = Z_T4_T5t-Z_T8_T9d+1; 

% lowerOffset = 15;
FH       = 18;  % widest slice of the formoral head
Z_L5_S1b = 57;  % bottom of L5-S1 intervertebral disk slice
%Z_L4_L5b = 69;   % bottom of L4 vertebra disk slice
Z_L4_L5  = 75;   % middle of L4-L5 intervertebral disk slice
Z_L3_L4  = 89;   % middle of L3-L4 intervertebral disk slice
UM       = 93;   % umbilicus (belly button)
Z_L2_L3  = 103;  % middle of L2-L3 intervertebral disk slice
Z_L1_L2  = 115;  % middle of L1-L2 intervertebral disk slice
%Z_L1_L2t    = 117; % top of L2 vertebra disk
lowerSlices = Z_L1_L2:-1:FH+2;
repeatEstimators = 0;

Nslices = length(upperSlices) + length(lowerSlices) + repeatEstimators;
slices = [upperSlices, lowerSlices];
levels(1:length(upperSlices))         = 1;
levels(length(upperSlices)+1:Nslices) = 2;
voxelVolumes(1:length(upperSlices))         = voxelVolume_upper;
voxelVolumes(length(upperSlices)+1:Nslices) = voxelVolume_lower;

% NOTE: belly button is from i = 130-122
bellybutton = 130:122;
bbx = 67;

% parameter settings
cutLineLength = 82;

backgroundThreshold   = 10;
vetoFactorF = 4.5;
vetoFactorW = 3.5;
foregroundThreshold1  = 59;
foregroundThreshold2  = 10;
foregroundThresholdF2 = 145;
foregroundThresholdV  = 103;
foregroundThresholdV2 = 230;

nbrThresholdsF = 7; %4; %7;
nbrThresholdsW = 7; %4; %9;
nbrThresholdsV = 7; %4; %7;

seSCAT = strel('rectangle',[7,17]);  
VATprelim      = 17;      % opening for preliminary VAT segmentation

waterXveto     = 73;
voidThreshold  = 5000;   % threshold of count of number of void pixels
lungErode      = 6;

aortaThreshold = 80;   % threshold of 8-bit water dominated signal data
aortaMinArea   = 70;
aortaDistance  = 50;
aortaRadius    = 9;

heartDilate    = 45;
notLungDilate  = 42; %21;
innerDilate    = 11; %13;
heartShift     = 25;
heartAreaThreshold = 400;

CATdilate      = 6; %3; %4; %6
CATmargin      = 6;
PAATdilate     = 6;

CATdilateCorrect = zeros(1,Nslices);
CATdilateCorrect(heartMax-23) = +1;
CATdilateCorrect(heartMax-24:heartMax-20) = +4;
CATdilateCorrect(heartMax-19) = +1;
CATdilateCorrect(heartMax-17:heartMax-13) = -1;
CATdilateCorrect(heartMax-1) = -1;
CATdilateCorrect(heartMax) = -1;
CATdilateCorrect(heartMax+1) = -2;
CATdilateCorrect(heartMax+2:heartMax+3) = -3;
CATdilateCorrect(heartMax+3:heartMax+4) = -2;
CATdilateCorrect(heartBottom-3:heartBottom-1) = +2;
CATdilateCorrect(heartBottom) = +3;

% load fat dominated signal data for upper scan
niiFup = load_nii(NIFTI_file_name_F_upper);
niiFlo = load_nii(NIFTI_file_name_F_lower);
% load water selected data for upper scan
niiWup = load_nii(NIFTI_file_name_W_upper);
niiWlo = load_nii(NIFTI_file_name_W_lower);

dims = size(niiFup.img);

% initializations
SCATvolume   = zeros(1,Nslices); % [mL]
VATvolume    = zeros(1,Nslices); % [mL]
organsVolume = zeros(1,Nslices); % [mL]
voidsVolume  = zeros(1,Nslices); % [mL]
lungVolume   = zeros(1,Nslices); % [mL]
heartVolume  = zeros(1,Nslices); % [mL]
aortaVolume  = zeros(1,Nslices); % [mL]
CATvolume    = zeros(1,Nslices); % [mL]
PAATvolume   = zeros(1,Nslices); % [mL]

aortaSeed = false(dims(2),dims(1),dims(3));
% <---- Add these here
% vertical section
%aortaSeed(,,)   = 1;
% curved section
%aortaSeed(,,)     = 1;




%% slice loop

for i = 1:Nslices
% must do for bottom of heart slices    
% for i = heartMax:heartBottom
% if i == heartBottom
%     showProgress = true;
%     batchShow    = ~showProgress;
% end
% problem slice:
%i = 16;
i 

slice = slices(i);    
level = levels(i);    
voxelVolume = voxelVolumes(i);    

%% make 8-bit slice image of fat

% orient image slice (posterior down)
if level == 1
    IFraw = fliplr(rot90(niiFup.img(:,:,slice)',2));
elseif level == 2
    IFraw = fliplr(rot90(niiFlo.img(:,:,slice)',2));
    IFraw = circshift(IFraw,[Y_shift,X_shift]);
end

% convert data to 8-bit gray-scale image
IFshow = uint8(double(IFraw)/double(max(max(IFraw)))*255);

if showProgress && ~batchShow
figure, imshow(IFshow);
end

IFrawVeto = (IFshow < vetoFactorF*backgroundThreshold);

%% make 8-bit slice image of water selected data

% orient image slice (posterior down)
if level == 1
    IWraw = fliplr(rot90(niiWup.img(:,:,slice)',2));
elseif level == 2
    IWraw = fliplr(rot90(niiWlo.img(:,:,slice)',2));
    IWraw = circshift(IWraw,[Y_shift,X_shift]);    
end

% convert data to 8-bit gray-scale image
IWshow = uint8(double(IWraw)/double(max(max(IWraw)))*255);

if showProgress && ~batchShow
    figure, imshow(IWshow);
end

IWrawVeto = (IWshow < vetoFactorW*backgroundThreshold);

%% correct fat signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
[IFshow1,IFcorrect] = CorrectInhomogeneity(IFraw,34); %100

if showProgress && ~batchShow
    figure, imshow(IFshow1);
end

%% segment the fat image foreground from background

% Find the peak and valley pixel bins of histogram
[peakPixelBinF,valleyPixelBinF] = FindHistogramExtrema(IFshow1.*uint8(1-IFrawVeto) + (IFshow.*uint8(IFrawVeto) + IFshow1.*uint8(IFrawVeto))/2,backgroundThreshold,foregroundThreshold1,(showProgress&&~batchShow));
peakPixelBinF
valleyPixelBinF

% rescale fat image using the histogram peak as the maximum  
peakIndicesF = find(IFshow1 == peakPixelBinF);
IFpeak = IFcorrect(peakIndicesF(1));
IFcorrect(IFcorrect>IFpeak) = IFpeak;
IFshow2 = uint8(double(IFcorrect)/double(IFpeak)*255);

if showProgress && ~batchShow
    figure, imshow(IFshow2);
end

% limit the valley pixel bin by Otsu's multithreshold method
valleyPixelBinF = AutorestrictThreshold(IFshow2,valleyPixelBinF,nbrThresholdsF);
valleyPixelBinF

% define fat foreground as being above histogram valley
seg_IFbin = (IFshow2 > valleyPixelBinF).*logical(1-IFrawVeto);

if showProgress && ~batchShow
figure, imshow(seg_IFbin,[]);
end

%% correct water signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
[IWshow2,IWcorrect] = CorrectInhomogeneity(IWraw,50); %100);

%% segment the image foreground from background within body cavity 

% mask corrected image with the internal (body cavity) mask
IWcorrect = double(IWshow2).*double(IWshow2>(IFshow2-4*backgroundThreshold));

% rescale image to the maximum intensity
IWmax = max(max(IWcorrect));
IWshow3 = uint8(double(IWcorrect)/double(IWmax)*255);

if showProgress && ~batchShow
    figure, imshow(IWshow3);
end

% save cross-section of saggital slice image
if i < Nslices - repeatEstimators
    Isaggital(:,Nslices-5-i+1) = IWshow3(:,166).*uint8(1-IWrawVeto(:,166));
end

% Find the peak and valley pixel bins of histogram
[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3.*uint8(1-IWrawVeto) + (IWshow.*uint8(IWrawVeto) + IWshow3.*uint8(IWrawVeto))/2,backgroundThreshold,foregroundThreshold1 - 2*backgroundThreshold,(showProgress&&~batchShow));
%[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(uint8((double(IWshow3).*((1+2.0*double(1-IWrawVeto))/3.0))),backgroundThreshold,foregroundThreshold1 - 2*backgroundThreshold,(showProgress&&~batchShow));
peakPixelBinW = peakPixelBinW - 2*backgroundThreshold
valleyPixelBinW = valleyPixelBinW

% rescale water image using the histogram peak as the maximum  
peakIndicesW = find(IWshow3 == peakPixelBinW);
IWpeak = IWcorrect(peakIndicesW(1));
IWcorrect(IWcorrect>IWpeak) = IWpeak;
IWshow4 = uint8(double(IWcorrect)/double(IWpeak)*255);

if showProgress && ~batchShow
    figure, imshow(IWshow4);
end

% limit the valley pixel bin by Otsu's multithreshold method
valleyPixelBinW = AutorestrictThreshold(IWshow4,valleyPixelBinW,nbrThresholdsW);
valleyPixelBinW

% define water foreground as being above histogram valley
seg_IWbin = (IWshow4 > valleyPixelBinW).*double(1-IWrawVeto);

if showProgress && ~batchShow
    figure, imshow(seg_IWbin,[]);
end

%% separation of SCAT from VAT

% define a cut line mask to enable filling of holes apart from body cavity
cutMask = logical(0*seg_IFbin);
cutMask(size(seg_IFbin,1)-cutLineLength:size(seg_IFbin,1),round(0.5*size(seg_IFbin,2))) = true;
cutSaveLine = bwmorph((cutMask.*seg_IFbin),'clean');

% link SCAT together for belly button
linkMask = logical(0*seg_IFbin);
if sum(bellybutton == i) == 1
    linkMask(bbx,round(0.5*(size(seg_IFbin,2)-cutLineLength)):round(0.5*(size(seg_IFbin,2)+cutLineLength))) = true;
end  

% preliminary body and perimeter mask
bodyMask = imfill(logical(1-FindLargestArea(IFrawVeto)+linkMask),'holes');
perimeterMask = logical(bwmorph(bodyMask,'dilate')-imerode(bodyMask,ones(7)));

filledFatBase = logical(imfill((~cutMask.*(bwmorph(logical(seg_IFbin.*logical(1-seg_IWbin)),'open')+linkMask+perimeterMask)),'holes')+cutSaveLine);

% preliminary separation of SCAT from VAT using morphological openning 
openedBWfat = imdilate(imerode(filledFatBase,ones(5)),ones(5));  

% construct mask of entire body and perimeter more properly
bodyMask = imfill(bwmorph(logical(openedBWfat+linkMask+seg_IWbin),'close'),'hole');
% construct mask of body perimeter
se1 = strel('rectangle',[11,11]);
perimeterMask = logical(bodyMask-imerode(bodyMask,se1));

% correct SCAT with inner cavity mask (and correct masks)
innerBase = imerode(imdilate(seg_IWbin.*(1-perimeterMask),ones(5)),ones(5));
% mask out the mammary tissue
if i < 22
    innerBase(1:waterXveto,1:round(0.4*size(seg_IWbin,2))) = 0;   %#ok<NASGU>
    innerBase(1:waterXveto,round(0.6*size(seg_IWbin,2)):end) = 0;   %#ok<NASGU>
    innerBase(1:waterXveto+10,1:round(0.2*size(seg_IWbin,2))) = 0;   %#ok<NASGU>
    innerBase(1:waterXveto+10,round(0.8*size(seg_IWbin,2)):end) = 0;   %#ok<NASGU>
end
% further corrections and constructing of the inner cavity mask
innerBase = FindLargestArea(innerBase,'multiple',0.06*voidThreshold);
innerBase = imerode(imfill(imdilate(innerBase,ones(7)),'holes'),ones(7));
innerMask = imdilate(imerode(innerBase,ones(5)),ones(5));

% define the inner perimeter mask of the SCAT from the inner cavity mask
SCATinnerMask = logical(imdilate(innerMask,se1)-innerMask);
SCATinnerMaskFill = bwmorph(imerode(imfill(imdilate(SCATinnerMask,ones(5)),'holes'),ones(5)),'close');
SCATinnerMask = logical(SCATinnerMaskFill - imerode(SCATinnerMaskFill,se1));
% correct inner mask
innerMask = logical(innerMask + imerode(SCATinnerMaskFill,se1));
%openedBWfat = openedBWfat.*logical(1-innerMask);

% apply cut line
filledBWfat = logical(imfill((~cutMask.*logical(openedBWfat+perimeterMask+SCATinnerMask)),'holes')+cutSaveLine);

% more thorough separation of SCAT from VAT using erosion preserving surface
erodedBWfat = logical(imerode(filledBWfat,se1)+perimeterMask);

% clean up SCAT base mask
SCATbase = bwmorph(erodedBWfat,'open');
SCATbase = FindLargestArea(SCATbase);

if showProgress && ~batchShow
    figure, imshow(SCATbase,[]);
end

%% segment SCAT (subcutaneaous adipose tissue)

% remove any small spurs and indentations from SCAT base mask
SCATclosed = bwmorph(bwmorph(SCATbase,'spur'),'close');

% dilate base mask to encompass "true" area 
SCATdilate = imdilate(SCATbase,seSCAT);

% construct SCAT mask as the union of fat foreground with dilated base
% and removing any small spurs or isolated pixels
SCATmask = bwmorph(bwmorph(seg_IFbin.*SCATdilate,'clean'),'spur');

% define resulting SCAT image
SCAT = uint8(SCATmask).*IFshow2;

if ~batchShow
    figure, imshow(SCAT,[]);
end

% preliminary segmentation of VAT 
VATmask = seg_IFbin.*(1-SCATmask);

% remove SCAT and VAT from water signal (organs)
seg_IWbin = seg_IWbin.*(1-SCATmask).*(1-VATmask);

%% segment water bearing organs

% clean and fill holes in water foreground for organ mask
organsMask = bwmorph(bwmorph(logical(seg_IWbin),'clean'),'fill');

%% correct the SCAT mask

% define skin mask
seSkin = strel('rectangle',[5,7]);
skinMask = imdilate(logical(bodyMask-bwmorph(bodyMask,'erode')),se1);

% define cut line for SCAT
cutSaveLine = bwmorph((cutMask.*seg_IFbin),'clean');

% use inner boundary mask to fill in SCAT properly
SCATfill = logical(imfill((logical(1-cutMask).*logical(SCATmask+SCATinnerMask)),'holes')+cutSaveLine).*seg_IFbin;
SCATcorrect = SCATfill.*logical(1-SCATmask);
SCATmask = logical(SCATmask + SCATcorrect);
SCATmask = FindLargestArea(SCATmask);
SCAT = uint8(SCATmask).*IFshow2;

% correct organs mask
organsMask = organsMask.*logical(1 - SCATcorrect);

%% segment internal mask (body cavity)
% fill in the SCAT mask for an isolated thoracic body mask
se3 = strel('square',3);
se7 = strel('square',7);
SCATbase2 = logical(imerode(imclose(SCATmask,se7),se3)+perimeterMask);
SCATfill = imfill(SCATbase2,'holes');

% define the fill base mask for further void and VAT segmentation
fillbase = SCATfill-SCATbase2;

% define the area outside the body cavity
externalMask = logical(1-innerMask);
if showProgress && ~batchShow
    figure, imshow(externalMask)
end
SCATplusExternalMask = bwmorph(bwmorph(seg_IFbin.*SCATdilate,'clean'),'spur')|externalMask;

%% segment voids (regions with low fat and water signal within the body cavity)

% define voids base as the regions not claimed by any existing mask
voidsBase = logical(1-seg_IWbin).*logical(1-organsMask).*logical(1-seg_IFbin).*logical(1-SCATplusExternalMask);

% define voids mask by opening its base
voidsMask = bwmorph(voidsBase,'open').*logical(1-organsMask).*logical(1-seg_IFbin).*logical(1-SCATplusExternalMask);
voids = uint8((peakPixelBinW/255.0)*(255-IWshow).*uint8(voidsMask));

%% segment VAT (visceral adipose tissue)

% define VAT mask
VATmask = VATmask.*innerMask.*logical(1-seg_IWbin);
VAT = uint8(VATmask).*IFshow2;
if ~batchShow
    figure, imshow(VAT,[]);
end

% define mask for non-fat areas
nonATmask = logical(1-SCATplusExternalMask-VATmask);
nonAT = uint8(nonATmask).*IFshow2;
if showProgress && ~batchShow
    figure, imshow(nonAT);
end

%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(voidsBase - voidsMask);
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF)); 
% add unassigne voxels to appropriate masks
VATcorrect = logical((voidsCorrect - organsCorrect).*innerMask);
SCATcorrect = logical((voidsCorrect - organsCorrect - VATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow.*uint8(organsMask);

if showProgress && ~batchShow
    figure, imshow(organs ,[]);
end

%% Lung segmentation
if i < diaphram + 1
%if    sum(sum(voidsMask)) > voidThreshold
    
% Combine equally weighted data sets
Iplus = sqrt((double(IFcorrect)/double(IFpeak)).^2 + (double(IWcorrect)/double(IWpeak)).^2);
IplusMax = max(max(Iplus));

% define "void" dominated image with in body cavity
IV = (IplusMax - Iplus).*double(innerMask);
IV(isnan(IV)) = 0;

%% correct void signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
[IVshow,~] = CorrectInhomogeneity(IV,150,'medfilt');

if showProgress && ~batchShow
    figure, imshow(IVshow);
end

% find peak and valley in histogram of the void image foreground
[peakPixelBinV,valleyPixelBinV] = FindHistogramExtrema(IVshow,backgroundThreshold,foregroundThresholdV,(showProgress&&~batchShow));
peakPixelBinV
valleyPixelBinV
% adjust valley back to original value
valleyPixelBinV = round(valleyPixelBinV*peakPixelBinV/255);
valleyPixelBinV

% restrict foreground threshold to be equal to or less than the next to
% highest automatically determined threshold
valleyPixelBinV = AutorestrictThreshold(IVshow,valleyPixelBinV,nbrThresholdsV);
valleyPixelBinV

% define void foreground base as being above histogram valley, use for lungs
lungBase = bwmorph((IVshow > valleyPixelBinV),'open');

% define lungMask
lungMask = bwmorph(FindLargestArea(lungBase,'multiple',0.1*voidThreshold),'close');

% correct existing masks
voidsCorrect = logical(uint8(voidsMask-lungMask));
voids(lungMask) = 0;
VATcorrect = ((VAT.*uint8(voidsCorrect))>=(organs.*uint8(voidsCorrect)))&voidsCorrect;
VATmask(VATcorrect) = 1;
VATmask = VATmask.*(1-lungMask);
VAT = uint8(VATmask).*IFshow2;
organsCorrect = logical(voidsCorrect - VATcorrect);
organsMask(organsCorrect) = 1;
organsMask = organsMask.*(1-lungMask);
organs = IWshow4.*uint8(organsMask);

%% aorta segmentation

% define the aorta base mask from seed
aortaBase = aortaSeed(:,:,slice); 

% grow aortaMask from seed by dilation
sea = strel('disk',aortaRadius);
aortaMask = imdilate(aortaBase,sea).*(IWshow>aortaThreshold);
aortaMask = imfill(aortaMask,'holes');

% choose largest area within field as aorta near the image center
[labeledBWaorta,nbr_labels_aorta] = bwlabel(aortaMask,4);
statsAorta = regionprops(labeledBWaorta,'area','centroid');
label_aorta_areas = cat(1, statsAorta.Area);
label_aorta_centroids = cat(1, statsAorta.Centroid);
k_aorta_max = max(label_aorta_areas);
aorta_maxLabel = find(label_aorta_areas == k_aorta_max);

% selection of largest area is restricted by certain conditions
if ~isempty(k_aorta_max)
    if k_aorta_max > aortaMinArea && ((label_aorta_centroids(aorta_maxLabel,1)-round(0.5*size(aortaMask,1))...
            + label_aorta_centroids(aorta_maxLabel,2)-round(0.5*size(aortaMask,2)))<aortaDistance);
        aortaMask = bwmorph((labeledBWaorta == aorta_maxLabel),'open');
    else
        aortaMask = 0*aortaMask;
    end
else
    aortaMask = 0*aortaMask;
end
aorta = IWshow4.*uint8(aortaMask);

% define aorta veto line for heart segmentation from aorta seed
aortaVeto = aortaBase;
sum(sum(aortaBase))
if sum(sum(aortaBase)) < 3
    aortaVeto = imdilate(aortaBase,ones(5,50));
else
    aortaBase = aortaSeed(:,:,slice-3); 
    aortaVeto = logical(1 + aortaMask);
end

%% heart segmentation

% restrict to only above observed bottom of heart
if i > heartBottom
    % do nothing
    heartMask = 0*organsMask;
else
    
% assume that heart is near the image center (shifted up and to the left)
heartBase = 0*organsMask;
heartBase(round(0.5*size(organsMask,1))-heartShift,round(0.5*size(organsMask,2))-heartShift) = 1;
seh = strel('disk',heartDilate);
heartBase = imdilate(heartBase,seh).*logical(organsMask+VATmask).*logical(1-aortaVeto);

% use lungs to frame the heart region
seNotlung = strel('disk',lungErode);
notLung = logical((1-lungMask).*fliplr(1-lungMask));
notLung = imerode(notLung.*heartBase,seNotlung);

% plot lung mask and centroid framed by lung
s = regionprops(notLung,'centroid');
centroids = cat(1, s.Centroid);

if showProgress && ~batchShow && ~isempty(centroids)
    figure, imshow(lungMask)
    hold on
    plot(centroids(:,1),centroids(:,2), 'b*')
    hold off
end

centroidMask = 0*notLung;
if ~isempty(centroids)
    centroidMask(round(centroids(:,2)),round(centroids(:,1))) = 1;
end

% remove side lobes that can be caused by diaphraim
lungVeto = bwmorph(imfill(imclose(lungMask,ones(35)),'holes')...
                   -imclose(lungMask,ones(35)),'dilate',6)...
                   .*logical(1-lungMask);
% unless it is the heart region ...              
if sum(sum(logical(lungVeto - lungVeto.*logical(1-centroidMask)))) == 1
    lungVeto = 0.*lungVeto;
end

% assume that the heart is framed by the lungs
notLung = FindLargestArea(notLung);
notLung = bwmorph(notLung.*logical(1-lungVeto),'open');

if showProgress && ~batchShow
    figure, imshow(notLung,[]);
end

% update heart base by dilation and masking
seHeart = strel('disk',notLungDilate);  %21
seInner = strel('disk',innerDilate);    %13
heartBase = imdilate(notLung,seHeart).*logical(1-lungMask)...
               .*logical(organsMask+VATmask)...
               .*logical(1-bwmorph(aortaMask,'dilate'))...
               .*logical(1-imdilate(SCATinnerMask,seInner))...
               .*logical(1-aortaVeto);
% get rid of spurs
heartBase = imerode(heartBase,ones(5));
           
% define the heart mask
heartMask = imdilate(FindLargestArea(heartBase),ones(5));
heartMask = bwmorph(heartMask.*logical(1-lungVeto),'open').*organsMask;
heartMask = FindLargestArea(heartMask,'multiple',heartAreaThreshold);

% mask lower part of heart with previous segmentation to restrict its area
if (i > heartMax) && (i < (heartBottom + 1))
    if i < (heartBottom-2)
      seHC = strel('disk',CATdilate-3);    
    else
      seHC = strel('disk',CATdilate-4);
    end
    heartCorrect = imerode(logical((heart3d(:,:,i-1)>0)+(CAT3d(:,:,i-1)>0)),seHC);
    heartCorrect = logical(heartCorrect + heartMask.*imdilate(VATmask,ones(7)) + heartMask.*imdilate(lungMask,ones(5)));
    heartCorrect = imdilate(imerode(heartCorrect,ones(5)),ones(5));
    heartMask = heartMask.*heartCorrect;
end

end

%% segment CAT ('X'-cardial adipose tissue near heart)

% define CAT via dilation
seCAT = strel('disk',CATdilate + CATdilateCorrect(i));
heartFatBase = bwmorph(bwmorph((imdilate(heartMask,seCAT).*VATmask)|heartMask,'spur'),'open');
CATbase = bwmorph(bwmorph(heartFatBase,'dilate'),'bridge');

% find largest fat area
CATmask = (FindLargestArea(CATbase).*VATmask); 

% correct heart mask with CAT mask
CATmaskSum = logical(sum(CATmask));
CATlen = length(CATmaskSum);
CATrightmost = CATlen;
for j = CATlen:-1:1
    if CATmaskSum(j) && (CATrightmost == CATlen)
        CATrightmost = j;
    end
end
CATveto = CATmask;
if CATrightmost < (CATlen-CATmargin)
    CATveto(:,1:CATrightmost+CATmargin) = 1;
else
    CATveto = ones(size(CATmask,1),size(CATmask,2));
end
heartMask = heartMask.*CATveto;
heartMask = FindLargestArea(heartMask,'multiple',heartAreaThreshold);

% define edge mask of void image to attempt to restrict to pericardial sack
CATsat = IFshow2.*uint8(bwmorph(CATmask,'erode'));
CATsat(CATsat==0) = 255;
CATedgeBase = bwmorph(128-CATsat,'clean');
CATedgeBase = CATedgeBase.*(1-imdilate(lungMask,ones(9)));
CATedgeBase = bwmorph(CATedgeBase,'clean');
seCEM = [[1,1,1,1,1];[0,1, 1, 1,0];[0,0,1,0,0];[0,0,0,0,0];[0,0,0,0,0]];
CATedgeMask = imdilate(CATedgeBase,seCEM);
CATmask = CATmask.*(1-CATedgeMask);

% fill in any holes in CAT and heart mask
CATheart = imfill(bwmorph(bwmorph(logical(CATmask+heartMask),'dilate',2),'erode',2),'holes');
heartMask = logical(heartMask + organsMask.*CATheart);
CATmask = logical(CATmask + VATmask.*CATheart);

%% segment PAAT (periaortic adipose tissue)

% define the PAAT base
sePAAT = strel('disk',PAATdilate);
aortaFatBase = logical(((imdilate(aortaMask,sePAAT).*VATmask)|aortaMask));
PAATbase = bwmorph(bwmorph(aortaFatBase,'dilate'),'bridge');
PAATmask = (FindLargestArea(PAATbase).*VATmask);

% correct PAATmask with aorta (require to not be above the aorta)
sePAATcorrect = ones(1+2*PAATdilate);
sePAATcorrect(1:PAATdilate,:) = 0;
PAATcorrect = imdilate(aortaMask,sePAATcorrect);
PAATmask = PAATmask.*PAATcorrect;

% correct CAT and VAT mask
CATmask = CATmask.*logical(1-PAATmask);
VATmask = VATmask.*logical(1-CATmask).*logical(1-PAATmask);
organsMask = organsMask.*logical(1-heartMask).*logical(1-aortaMask);

%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(bwmorph(bodyMask,'erode',3) - logical(SCATmask + VATmask + voidsMask));
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF)); 
% add unassigne voxels to appropriate masks
VATcorrect = logical((voidsCorrect - organsCorrect).*innerMask);
SCATcorrect = logical((voidsCorrect - organsCorrect - VATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow.*uint8(organsMask);

% adipose tissue array
lung     = IVshow.*uint8(lungMask);
voids    = lung;
ATarray  = [SCAT,VAT,organs,voids];

tissues(:,:,1) = organs + SCAT;
tissues(:,:,2) = VAT + SCAT;
tissues(:,:,3) = voids;

if showProgress && ~batchShow
    figure, imshow(tissues);
end

SCATarea = sum(sum(SCATmask))
VATarea  = sum(sum(VATmask))
organArea  = sum(sum(organsMask))
voidArea = sum(sum(voidsMask))

SCATfuzzyArea   = sum(sum(SCAT))/255.
VATfuzzyArea    = sum(sum(VAT))/255.
organFuzzyArea  = sum(sum(organs))/255.
voidFuzzyArea   = sum(sum(voids))/255.

% thoracic tisse array
heart    = IWshow4.*uint8(heartMask);
CAT = IFshow2.*uint8(CATmask);
PAAT = IFshow2.*uint8(PAATmask);

% lung-heart-fat array
TTarray = [lung,heart,CAT];
%if showProgress && ~batchShow
%figure, montage(TTarray);
%end

TTtissues(:,:,1) = heart + aorta + CAT;
TTtissues(:,:,2) = VAT;
TTtissues(:,:,3) = lung + aorta + PAAT + CAT;

figure, imshow(TTtissues);

heartTissues(:,:,1) = heart + CAT +aorta;
heartTissues(:,:,2) = CAT + PAAT;
heartTissues(:,:,3) = aorta + PAAT + CAT;

if ~batchShow
figure, imshow(heartTissues);
end

SCATvolume(i)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
VATvolume(i)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
organsVolume(i)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
voidsVolume(i) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]
lungVolume(i)  = double(sum(sum(lungMask)))*voxelVolume/1000.0;  % [mL]
heartVolume(i) = double(sum(sum(heartMask)))*voxelVolume/1000.0; % [mL]
aortaVolume(i) = double(sum(sum(aortaMask)))*voxelVolume/1000.0; % [mL]
CATvolume(i)   = double(sum(sum(CATmask)))*voxelVolume/1000.0;   % [mL]
PAATvolume(i)  = double(sum(sum(PAATmask)))*voxelVolume/1000.0;  % [mL]

SCAT3d(:,:,i)   = SCAT;
VAT3d(:,:,i)    = VAT;
organs3d(:,:,i) = organs;
voids3d(:,:,i)  = voids;
lung3d(:,:,i)   = lung;
heart3d(:,:,i)  = heart;
aorta3d(:,:,i)  = aorta;
CAT3d(:,:,i)    = CAT;
PAAT3d(:,:,i)   = PAAT;

else

%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(bwmorph(bodyMask,'erode',3) - logical(SCATmask + VATmask + voidsMask));
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF)); 
% add unassigne voxels to appropriate masks
VATcorrect = logical((voidsCorrect - organsCorrect).*innerMask);
SCATcorrect = logical((voidsCorrect - organsCorrect - VATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow.*uint8(organsMask);    
    
% adipose tissue array
ATarray = [SCAT,VAT,organs,voids];

tissues(:,:,1) = organs + SCAT;
tissues(:,:,2) = VAT + SCAT;
tissues(:,:,3) = voids;

figure, imshow(tissues);

if i == (length(upperSlices) + 1)

SCATvolume(i-1)  = 0.5*(SCATvolume(i-1) + double(sum(sum(SCATmask)))*voxelVolume/1000.0);  % [mL]
VATvolume(i-1)   = 0.5*(VATvolume(i-1)  + double(sum(sum(VATmask)))*voxelVolume/1000.0);   % [mL]
organsVolume(i-1)= 0.5*(organsVolume(i-1)+double(sum(sum(organsMask)))*voxelVolume/1000.0);% [mL]
voidsVolume(i-1) = 0.5*(voidsVolume(i-1) +double(sum(sum(voidsMask)))*voxelVolume/1000.0); % [mL]

SCAT3d(:,:,i-1)  = 0.5*(SCAT3d(:,:,i-1)+SCAT);
VAT3d(:,:,i-1)   = 0.5*(VAT3d(:,:,i-1)+VAT);
organs3d(:,:,i-1)= 0.5*(organs3d(:,:,i-1)+organs);
voids3d(:,:,i-1) = 0.5*(voids3d(:,:,i-1)+voids);

elseif i > (length(upperSlices) + 1)
    
SCATvolume(i-1)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
VATvolume(i-1)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
organsVolume(i-1)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
voidsVolume(i-1) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]

SCAT3d(:,:,i-1)  = SCAT;
VAT3d(:,:,i-1)   = VAT;
organs3d(:,:,i-1)= organs;
voids3d(:,:,i-1) = voids;

else

SCATvolume(i)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
VATvolume(i)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
organsVolume(i)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
voidsVolume(i) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]

SCAT3d(:,:,i)  = SCAT;
VAT3d(:,:,i)   = VAT;
organs3d(:,:,i)= organs;
voids3d(:,:,i) = voids;    
    
end

end

end

% properly orient resulting saggital image
Isaggital = rot90(Isaggital',2); 

% save results
save('subject01_i.mat')