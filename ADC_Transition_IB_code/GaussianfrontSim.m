
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for generating model "cell density" maps with gaussian "infiltration fronts" and corresponding
% diffusion maps
%
% Dr Antoine Vallatos / Dr Gerry Thompson
% University of Edinburgh
%
% Version 1.6
% Tuesday 14th Oct 2023

% MATLAB Toolboxes required: 
%Image Prossessing

%External code required:
%NIFTI analysis: Jimmy Shen (2023). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. 
%Multi Frame Images: Felipe G. Nievinski (2023). subtightplot (https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot), MATLAB Central File Exchange. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
set(0,'defaultfigurecolor',[1 1 1])

%%  Define input output folders
% Run and select the mouse folder (no need to go in the folder)
uiwait(helpdlg('Select input folder',...
    ' data analysis'));
inputfolder = uigetdir
outputfolder = strcat(inputfolder, sprintf('\\Front_Analysis\')) %USER ACTION Define folder
mkdir(outputfolder)

%% Generate 2D Gaussian distribution
% At the resolution of MRI
close all;
size=176; %preclinical MRI image resolution
g=zeros(size,size); %2D filter matrix
l=zeros(size,size);
sigma=16; %standard deviation of the gaussian distribution (~invasiveness)- values 8 16 or 32
normdens=0.1 %percentage of normal cells
Vplateaumax=0.9
Vplateaumin=0.01

%Gaussian & lin filter
for i=-(size-1)/2:(size-1)/2
    for j=-(size-1)/2:(size-1)/2
        x0=(size+1)/2; %center
        y0=(size+1)/2; %center
        x=i+x0; %row
        y=j+y0; %col
        g(y,x)=exp(-((x-x0)^2+(y-y0)^2)/2/sigma/sigma);
        l(y,x)= sqrt((x-x0)^2+(y-y0)^2);
    end
end

%Plot position map
figure('Name','Position','NumberTitle','off');
set(gcf, 'Position', [110 500 400 400],'color','w')
imagesc(l);
colormap(gray);
saveas(gcf,strcat(outputfolder, sprintf('\\lin')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\lin')), 'pdf');
figure; plot(l(:,size/2));
set(gcf, 'Position', [110 10 400 400],'color','w')

g=g+normdens
gpl=g;
gpl(gpl>=Vplateaumax)=Vplateaumax;
gpl(gpl<=Vplateaumin)=Vplateaumin;

%Gaussian density distribution
figure('Name','Density','NumberTitle','off');
set(gcf, 'Position', [610 500 400 400],'color','w')
imagesc(gpl);
colormap(gray);
saveas(gcf,strcat(outputfolder, sprintf('\\dens')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\dens')), 'pdf');

% Density with plateau at 10% higher values 
figure('Name','Density plateau','NumberTitle','off');
set(gcf, 'Position', [1010 500 400 400],'color','w')
imagesc(gpl);
colormap(gray);
saveas(gcf,strcat(outputfolder, sprintf('\\denspl')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\denspl')), 'pdf');

% Profile vector
figure; plot(g(:,size/2));
set(gcf, 'Position', [1410 500 400 400],'color','w')
hold
plot(gpl(:,size/2));
xlim([0 size])
ylim([0 1.2])


%%Generate diffusion maps

% Normalize gaussian filter 
sum1=sum(g);
sum2=sum(sum1);
g=g/sum2;
sum1=sum(gpl);
sum2=sum(sum1);
gpl=gpl/sum2;


%Calculate Diffusion values
diff=1./sqrt(g);
diffpl=1./sqrt(gpl);
diffrg=diff;
diffrgpl=diffpl;
diffrg=diffrg*(0.85e-3/max(max(diffrg)));
diffrgpl=diffrgpl*(0.85e-3/max(max(diffrgpl))); % ADC plateau 0.9x10-9 mm2/s

% plot gaussian density diffusion distribution
figure('Name','Diffusion reg','NumberTitle','off');
set(gcf, 'Position', [610 10 400 400],'color','w')
imagesc(diffrg);
colormap(gray);colorbar
saveas(gcf,strcat(outputfolder, sprintf('\\diffrg')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\diffrg')), 'pdf');

% plot diffusion distribution for density with plateau for 10% higher density (~ "tumour core")
figure('Name','Diffusion plat reg','NumberTitle','off');
set(gcf, 'Position', [1010 10 400 400],'color','w')
imagesc(diffrgpl);
colormap(gray);colorbar;
saveas(gcf,strcat(outputfolder, sprintf('\\diffrgplat')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\diffrgplat')), 'pdf');

% plot profiles
figure; plot(diffrg(:,size/2));
set(gcf, 'Position', [1410 10 400 400],'color','w')
hold
plot(diffrgpl(:,size/2));

%% Density ROIS
Thresh_coeff=(max(g(:))-min(g(:)))/5
thresholds=[5 4 3 2 1]*Thresh_coeff 

figure('Name','Thresholds','NumberTitle','off');
set(gcf, 'Position', [0 0 1000 200],'color','w')
for i=1:5
    
    gth{i}=g;
    gth{i}(gth{i}<=thresholds(i))=0;
    gth{i}(gth{i}>thresholds(i))=1;
    subtightplot(1,5,i,0.01);
    imagesc(gth{i});axis off
    colormap(gray);
    nb(i)=sum(sum(gth{i}))
end
saveas(gcf,strcat(outputfolder, sprintf('\\gaussthresh')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\gaussthresh')), 'pdf');

% Pixel count per ROI
figure('Name','Rois','NumberTitle','off');
set(gcf, 'Position', [110 100 1000 200],'color','w')
for i=5:-1:1
    if i==1
        groi{i}=gth{i};
        subtightplot(1,5,i,0.01);
        imagesc(groi{i});colormap(gray);axis off
        
    else
        
        groi{i}=gth{i}-gth{i-1}
        subtightplot(1,5,i,0.01);
        imagesc(groi{i});colormap(gray);axis off
        
    end
    
    countg(i)=sum(sum(groi{i}(groi{i}~=0)));
end
saveas(gcf,strcat(outputfolder, sprintf('\\gaussroi')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\gaussroi')), 'pdf');

figure('Name','Counts','NumberTitle','off');plot(countg);
saveas(gcf,strcat(outputfolder, sprintf('\\countroielems')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\countroielems')), 'pdf');

%% Density per ROI
close all;

%Plot values for density map based masks
figure('Name','DensMasks','NumberTitle','off');
set(gcf, 'Position', [110 500 1000 200],'color','w')
for i=1:5
    gmask{i}=groi{i}.*g
    subtightplot(1,5,i,0.01);
    imagesc(gmask{i});colormap(gray);axis off
    meang(i)=mean2(mean2(gmask{i}(gmask{i}~=0)));
end
saveas(gcf,strcat(outputfolder, sprintf('\\gaussmask')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\gaussmask')), 'pdf');

%Plot density / mask
figure('Name','DensMean','NumberTitle','off');plot(meang);
set(gcf, 'Position', [1210 500 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'pdf');

%Plot values for density map (with plateau) based masks
figure('Name','DensMasksplateau','NumberTitle','off');
set(gcf, 'Position', [110 100 1000 200],'color','w')
for i=1:5
    gmaskpl{i}=groi{i}.*gpl
    subtightplot(1,5,i,0.01);
    imagesc(gmaskpl{i});colormap(gray);axis off
    meangpl(i)=mean2(mean2(gmaskpl{i}(gmaskpl{i}~=0)));
end
saveas(gcf,strcat(outputfolder, sprintf('\\gaussmaskplat')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\gaussmaskplat')), 'pdf');

%Plot density / mask (w plateau)
figure('Name','DensMean','NumberTitle','off');plot(meangpl);
set(gcf, 'Position', [1210 100 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffplatrois')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffplatrois')), 'pdf');

%% Diffusion per ROI
close all;

%Plot diffusion values for density map based masks
figure('Name','diffMasks','NumberTitle','off');
set(gcf, 'Position', [110 500 1000 200],'color','w')
for i=1:5
    diffmask{i}=groi{i}.*diffrg
    subtightplot(1,5,i,0.01);
    imagesc(diffmask{i});colormap(gray);axis off
    meandiff(i)=mean2(mean2(diffmask{i}(diffmask{i}~=0)));
end
saveas(gcf,strcat(outputfolder, sprintf('\\diffmask')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\diffmask')), 'pdf');

%Plot diffusion / mask
figure('Name','DiffMean','NumberTitle','off');plot(meandiff);
set(gcf, 'Position', [1210 500 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'pdf');


%Plot diffusion values for density map (with plateau) based masks
figure('Name','diffMasksplateau','NumberTitle','off');
set(gcf, 'Position', [110 100 1000 200],'color','w')
for i=1:5
    diffmaskpl{i}=groi{i}.*diffrgpl
    subtightplot(1,5,i,0.01);
    imagesc(diffmaskpl{i});colormap(gray);axis off
    meandiffpl(i)=mean2(mean2(diffmaskpl{i}(diffmaskpl{i}~=0)));
end
saveas(gcf,strcat(outputfolder, sprintf('\\diffmask')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\diffmask')), 'pdf');

%Plot diffusion / mask(w plateau)
figure('Name','DiffMeanpl','NumberTitle','off');plot(meandiffpl);
set(gcf, 'Position', [1210 100 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffplrois')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffplrois')), 'pdf');



%% Diffusion vs Density plots
%Plot diffusion / density
figure('Name','DiffvsDens','NumberTitle','off');plot(meang,meandiff);
set(gcf, 'Position', [1510 500 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffrois')), 'pdf');
%Plot diffusion / density (w plateau)
figure('Name','DiffvsDensplat','NumberTitle','off');plot(meangpl,meandiffpl);
set(gcf, 'Position', [1510 100 300 300],'color','w')
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffroispl')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\meandiffroispl')), 'pdf');

%% Export for Vectorial analysis code - 
%% !!!!!!!!!!!!!!!!Go to a reference mouse folder!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%e.g. M05
% ROIs attributed to ADC, T2W and HLA regions with ROIADC<ROIT2W<ROIHLA
ADCroidata=groi{1}(:,:)+groi{2}(:,:);
T2Wroidata=groi{1}(:,:)+groi{2}(:,:)+groi{3}(:,:);
HLAroidata=groi{1}(:,:)+groi{2}(:,:)+groi{3}(:,:)+groi{4}(:,:);

% Replace ADC, T2W and HLA data in importered mouse dcm data (e.g. M10)
ADCROI=load_nii('ROIAVADC.nii')
ADCROI.img=ADCroidata
T2WROI=load_nii('ROIAVT2W.nii')
T2WROI.img=T2Wroidata
HLAROI=load_nii('ROIAVHLA.nii')
HLAROI.img=HLAroidata
M08=load_untouch_nii('M05wk12.nii.gz')%Load reference animal

%% Create a pseudo-mouse dataset
%% !!!!!!!!!!!!Go to the folder of the newly created "mouse"!!!!!!!!!!!!!!!!!!!!!!!!

% Save ROIs and replace diffusion data by diffusion plateau data
save_nii(ADCROI,'ROIAVADC.nii');
save_nii(T2WROI,'ROIAVT2W.nii');
save_nii(HLAROI,'ROIAVHLA.nii');

save("ADCroidata.mat","ADCroidata");
save("T2Wroidata.mat","T2Wroidata");
save("HLAroidata.mat","HLAroidata");

volume_image=zeros(176,176,4);
for i=1:2
    volume_image(:,:,i)=diffrgpl;
end
for i=3:4
    volume_image(:,:,i)=gpl;
end


save('matlabM11wk12.mat', 'volume_image');
% Save pseudo-mouse based on reference mouse dataset
M11=M01;
M11.img=volume_image;
save_untouch_nii(M11,'M11wk12.nii.gz'); %Save image as a model preclinical dataset (e.g. M11 week 12)

% Save Workspace
save(strcat(outputfolder, sprintf('\\Workspace.mat'))); % Type in name of file.

%% Front Profile plot

%Density profile

figure('Name','Density profile','NumberTitle','off');
plot(l,gpl,'k');axis([0 130 0 1.1*max(gpl(:))])
saveas(gcf,strcat(outputfolder, sprintf('\\densslope')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\densslope')), 'pdf');

%Diffusion profile
figure('Name','Diffusion profile','NumberTitle','off');
plot(l,diffrgpl, 'r'); axis([0 130 0.9*min(diffrgpl(:)) 1.1*max(diffrgpl(:))]);
saveas(gcf,strcat(outputfolder, sprintf('\\diffslope')), 'tif');
saveas(gcf,strcat(outputfolder, sprintf('\\diffslope')), 'pdf');

