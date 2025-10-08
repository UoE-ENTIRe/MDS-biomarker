
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for exploring 2D tumour boundary profiles in multispectral,
% co-registered mouse data from Bruker scanners

% Copyright
% Dr Gerard Thompson & Dr. Antoine Vallatos
% University of Edinburgh
% gerard.thompson@ed.ac.uk
%Version 2.1
% Tuesday 15th Oct 2023

% MATLAB Toolboxes required: 
%Image Prossessing, Fuzzy Logic

%External code required:
%NIFTI analysis: Jimmy Shen (2023). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. 
%Errorbars: Jean-Yves Tinevez (2023). errorbarxy (https://www.mathworks.com/matlabcentral/fileexchange/19002-errorbarxy), MATLAB Central File Exchange. 
%Overlay Images: Matt Smith (2023). imoverlay (https://www.mathworks.com/matlabcentral/fileexchange/42904-imoverlay), MATLAB Central File Exchange.
%Multi Frame Images: Felipe G. Nievinski (2023). subtightplot (https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot), MATLAB Central File Exchange. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
set(0,'defaultfigurecolor',[1 1 1])

% %%  Define input output folders
uiwait(helpdlg('Select input folder',...
    ' data analysis'));
inputfolder = uigetdir
outputfolder = strcat(inputfolder, sprintf('\\Vector_Analysis_Results_Sim\')) %USER ACTION Define folder!!
mkdir(outputfolder)

%% Initiate parameters
% Switch from function to script, make 1 or 0
itsascript = 1;
%initiate variables
if itsascript
    mouselist=[];
    k=[];
    sampletype=[];
    prof=[];
    f=[];
    interptype=[];
    moatsize=[];
end

%Comment to analyse simulated data
% mouselist =[1,3:10,14]; %mouse datasets to include in the analysis - no data for mouse 2 / mouse 14 is mouse 9 with HLA eplaced by HnE

% Uncomment to analyse Simulated data
mouselist =[11:13]; %model datasets to include in the analysis

%populating the data matrix (data_proc)
saveoutslices = 0;
saveout = 1;
loadall = 1; %If this condition is true, it means that there is a requirement to load all data.
loaddata_proc = 0; %If this condition is true, it means that the data_proc variable has not been loaded from any external source.
savedata_proc = 1; %If this condition is true, it means that some data has been saved in the data_proc variable.
pdims = [0.11 0.11 1.5]; %depends on imaging but already harmonised on this dataset 

% This is set to populate a multispectral sampling result matrix
numconts = 4;
plotshow = 1;
scatshow = 1;
clusshow = 1;
boxpshow = 1;

% Data provided here as a 3D  stack of a single slice
% 1.T2w 2.ADC 3.HLA 4.H&E. (M03, M07, M08 no H&E).

% For figures etc...
t2wnum = 1;
adcnum = 2;
hlanum = 3;
hnenum = 4;

labelvec = ['T2w','ADC','HLA','H&E'];

%%%% Load all datasets file
load('data_proc.mat');
% Mouse data: 1-10
% Modelled data 11-13
%Mouse 14 = Mouse 9 with HLA swap for H&E - see supplemental info
data_proc{14}=data_proc{9}
data_proc{14}(:,:,3)=data_proc{14}(:,:,4); 


%% Define Parameters
% ROIfile, mouselist, moat.
% Load the tumour ROIs
MASKS = cell(1);
SEGS = cell(1);
BW1 = cell(1);
BW1wide = cell(1);
BW2 = cell(1);

% roifiles to pick in eack mouse folder
roifile = ('./ROIAVADC.nii');
HLAroifile = ('./ROIAVHLA.nii');
T2Wroifile = ('./ROIAVT2W.nii');

% mouselist = input('Please input the mouse list you want to process(default is 1:10) ->: ')

if isempty(mouselist)
    % mouselist = (1:10);
    mouselist = (1:2);
end

% moatsize = input('Please input the moat size (number of dilations) you want (default is 4) ->: ');
if isempty(moatsize)
    moatsize = (8);
end

%% Surface normal generation (pseudo surface with padding)

% k = input('Please input the kernel size you want (default is 5) ->: ');
if isempty(k)
    k = 5;
end
qq = (k+1)/2;

% sampletype = input('Please state sampling method: 1 = boundary outwards | 2 = across boundary from inside to out (default 2) ->: ');
if isempty(sampletype)
    sampletype = 2;
end

% prof = input('Please input the length of the sampling profile (default = 6 voxels) ->: ');
if isempty(prof)
    prof = 6;
end
1
% f    = input('Please input the sampling frequency (default = 2 per voxel) ->: ');
if isempty(f)
    f = 2;
end
step = prof*f;
error = 1;

while error ~= 0
    % interpolation is provided since the sampling vectors must be oversampled to account for non-orthogonal trajectory and oblique traversing of voxels...
    % interptype = input('Please state interpolation method: 1 = nearest | 2 = linear | 3 = spline | 4 = cubic (default linear) ->: ');
    if isempty(interptype)
        interptype = 2;
    end
    switch interptype
        case 1
            interpstring = '*nearest';error = 0;
        case 2
            interpstring = '*linear';error = 0;
        case 3
            interpstring = '*spline';error = 0;
        case 4
            interpstring = '*cubic';error = 0;
        otherwise
            disp('Please choose again');error = 1;
    end
end

%% Loop through datasets

for xxx = mouselist
    dirname = sprintf('M%02dwk12',xxx);
    cd(dirname);
    
    %getting masks for mice and modelled data
    if (xxx~=11)&&(xxx~=12)&&(xxx~=13)
   
    ROI = load_untouch_nii('ROIAVADC.nii')
    HLAROI=load_untouch_nii('ROIAVHLA.nii')
    T2WROI = load_untouch_nii('ROIAVT2W.nii')
    
    ROI.img=flip(ROI.img,1); 
    HLAROI.img=flip(HLAROI.img,1); 
    T2WROI.img=flip(T2WROI.img,1); 

    MASKS{xxx} = ROI.img;
    HLAMASKS{xxx} = HLAROI.img;
    T2WMASKS{xxx} = T2WROI.img;
    else 
        
    load("ADCroidata.mat");
    load("HLAroidata.mat");
    load("T2Wroidata.mat");
        
    MASKS{xxx} = ADCroidata;
    HLAMASKS{xxx} = HLAroidata;
    T2WMASKS{xxx} =  T2Wroidata; 
    end
    
    SEGS{xxx} = int8(imbinarize(data_proc{xxx}(:,:,2))); % make a crude "brain mask" from ADC slice 9
    brain = make_nii(SEGS{xxx}, pdims); % save it out as NIfTI
    save_nii(brain,'BrainMask.nii.gz');
    
   % ROI boundary generation
    BW1{xxx} = edge(MASKS{xxx},'sobel');
    BW2{xxx} = edge(MASKS{xxx},'canny');
    
    % Dilate it for assessing the boundary moat/band
    BW1wide{xxx} = zeros([size(BW1{xxx}) moatsize]);
    for aa = 1:moatsize
        se = strel('square',aa);
        BW1wide{xxx}(:,:,aa) = imdilate(BW1{xxx},se);
    end
    
    cd ..
end


%% Histology normalisation
% 
% % HLA norm
% close all;
% for xxx = mouselist
%     C=squeeze(data_proc{xxx}(:,:,3));
%     C(SEGS{xxx}==0)=0;
%     Ccol = C(:);
%     pctl = 5:5:100;
%     for k1 = 1:length(pctl)
%         pctval{xxx}(k1,:) = [pctl(k1) prctile(Ccol(Ccol~=0),pctl(k1))];
%     end
%     
%     C(C<pctval{xxx}(2,2))=pctval{xxx}(2,2);
%     C(C>pctval{xxx}(19,2))=pctval{xxx}(19,2);
%     C=(C-min(min(C)))./(max(max(C))-min(min(C)));
%     data_proc{xxx}(:,:,3)=C;
%     
% end
% 
% % HnE norm
% for xxx = mouselist
%     C2=squeeze(data_proc{xxx}(:,:,4));
%     C2(SEGS{xxx}==0)=0;  
%     C2col = C2(:);
%     pctl = 5:5:100;
%     for k1 = 1:length(pctl)
%         pctval2{xxx}(k1,:) = [pctl(k1) prctile(C2col(C2col~=0),pctl(k1))];
%     end  
%     C2(C2<pctval2{xxx}(1,2))=pctval2{xxx}(1,2);
%     C2(C2>pctval2{xxx}(19,2))=pctval2{xxx}(19,2);
%     C2=(C2-min(min(C2)))./(max(max(C2))-min(min(C2))); 
%     data_proc{xxx}(:,:,4)=C2;
%     
% end

%% Profiles

tu = 1;

for xxx = mouselist
    clear ft fu fv fvcorr voxellist
    dirname = sprintf('M%02dwk12',xxx);
    cd(dirname);
    VOIedge = BW1{xxx};
    VOI = MASKS{xxx};
    ol = 21; %this is the padding to avoid any edge issues as unlikely kernel >7 in size
    % Make blanks to receive the block (filled) and ribbon (unfilled) VOIs
    VOIup = zeros(size(VOI,1),size(VOI,2),ol);
    VOIfl = zeros(size(VOI,1),size(VOI,2),ol);
    for u = 1:ol
        VOIup(:,:,u) = VOIedge;
        VOIfl(:,:,u) = VOI;
    end
    Ribn = make_nii(VOIup, pdims);
    save_nii(Ribn,'RIBBON.nii.gz');
    Blok = make_nii(VOIfl, pdims);
    save_nii(Blok,'BLOCK.nii.gz');
    [voxellist(:,1),voxellist(:,2)] = ind2sub(size(VOIedge),find(VOIedge>0));
    voxellist(:,3) = (ol+1)/2;
    % determine a centroid/centre of mass for fixing profile direction (caution
    % if very convoluted surface!
    voicentroidxx = (min(voxellist(:,1))+(max(voxellist(:,1)-min(voxellist(:,1))))/2);
    voicentroidyy = (min(voxellist(:,2))+(max(voxellist(:,2)-min(voxellist(:,2))))/2);
    voicentroidzz = (ol)/2;
    voicentroidx  = round(min(voxellist(:,1))+(max(voxellist(:,1)-min(voxellist(:,1))))/2);
    voicentroidy  = round(min(voxellist(:,2))+(max(voxellist(:,2)-min(voxellist(:,2))))/2);
    voicentroidz  = (ol+1)/2;
    centroid      = [voicentroidxx voicentroidyy voicentroidzz];
    ct            = [voicentroidx voicentroidy voicentroidz];
    voxellisttemp = voxellist;
    MASK = VOIup;
    SE1 = zeros(size(MASK,1),size(MASK,2),size(MASK,3),3,'single');
    SE2 = zeros(size(MASK,1),size(MASK,2),size(MASK,3),3,'single');
    SE3 = zeros(size(MASK,1),size(MASK,2),size(MASK,3),3,'single');
    % Make the blank kernel and surface voxel vector lists...
    p = zeros(k,k,k);
    q = zeros(k,k,k);
    s = [k,k,k];
    ft = zeros([length(voxellist) 3]);
    fu = zeros([length(voxellist) 3]);
    fv = zeros([length(voxellist) 3]);
    
    h = waitbar(0,sprintf('Surface Tangent Plane Fitting for Mouse %02d',xxx));
    for i = 1:length(voxellist)
        waitbar(i/length(voxellist));
        % populate the kernel
        for d = 1:k
            for e = 1:k
                for f = 1:k
                    p(d,e,f) = VOIup((voxellist(i,1)+d-qq),(voxellist(i,2)+e-qq),(voxellist(i,3)+f-qq));
                    q(d,e,f) = VOIfl((voxellist(i,1)+d-qq),(voxellist(i,2)+e-qq),(voxellist(i,3)+f-qq));
                end
            end
        end
        % calculate the eigenvectors relating to the 'surface'
        [Q,R,S] = ind2sub(s,find(p>0));
        P = [Q R S];
        CM = cov(P);
   
        [V,D] = eig(CM);
        
        ft(i,:) = V(:,1);
        fu(i,:) = V(:,2);
        fv(i,:) = V(:,3);
        
        
        % Ones or zeros sampling
        trm = 3;
        gp = zeros(trm,1);
        gm = zeros(trm,1);
        for ss = 1:trm
            tempVp = voxellist(i,:) + ss*ft(i,:);
            tempVm = voxellist(i,:) - ss*ft(i,:);
            gp(ss) = VOIfl(ceil(tempVp(1)),ceil(tempVp(2)));
            gm(ss) = VOIfl(ceil(tempVm(1)),ceil(tempVm(2)));
        end
        
        % in theory, more likely to sample all ones going centripetally/inwards
        if sum(gp) > sum(gm)
            ft(i,:) = -ft(i,:);
        end
        
    end
    
    close(h)
    
    voxellist = voxellisttemp;
    
    if plotshow == 1
 
    end
    
    % Generation of profiles and derived metrics
    % Set out some matrices to produce images for comparison
    traces = zeros(size(VOI));
    pointindex = zeros(size(VOI));
    lines = zeros(size(VOI));
    pvalscompnc = zeros(size(VOI));
    pvalscompc = zeros(size(VOI));
    corrcoef = zeros(size(VOI));
    ADCtcoef = zeros(size(VOI));
    HLAtcoef = zeros(size(VOI));
    HnEtcoef = zeros(size(VOI));
    

    % The spatially corrected vector is simply therefore the unit vector...
    % fvcorr = ft; NOTE this is fv in the 3D model 
    fvcorr = ft(:,1:2);
    
    % get the voxellist back to 2d...
    voxellist = voxellist(:,1:2);
    
    % initialise the sampling result matrix
    sample_mat = NaN(length(voxellist),step+1,numconts);
    voxelvalue = [0 0];
    
    % The sampling loops for starting at the edge
    if sampletype == 1
        fprintf('Measuring the voxel values (method %d) using %s interpolation for mouse %02d...\n',sampletype,interpstring(2:end),xxx);
        w = waitbar(0,sprintf('Measuring the voxel values (method %d) using %s interpolation for mouse %02d...',sampletype,interpstring(2:end),xxx));
        for i = 1:length(voxellist)
            waitbar(i/length(voxellist));
            pointindex(voxellist(i,1),voxellist(i,2)) = i;
            for q = 0:step
                p = q+1;
                voxelvalue = (voxellist(i,:)+((q/f).*fvcorr(i,:)));
                x = voxelvalue(1);
                y = voxelvalue(2);
                % avoid going out of bounds or sampling excluded tissue
                % classes
                if (1 < x) && (x < size(MASK,1)-1) && (1 < y) && (y < size(MASK,2)-1)
                    if SEGS{xxx}(ceil(x),ceil(y)) > 0
                        for b = 1:numconts
                            sample_mat(i,p,b)  = interpn(data_proc{xxx}(:,:,b),x,y,interpstring);
                        end
                        traces(ceil(x),ceil(y)) = traces(ceil(x),ceil(y)) + 1;
                        % Check direction of travel
                        lines(ceil(x),ceil(y)) = p;
                    else
                        break
                    end
                else
                    break
                end
            end
        end
        % The sampling loops for sampling across the edge
    elseif sampletype == 2
        fprintf('Measuring the voxel values (method %d) using %s interpolation for mouse %02d...\n',sampletype,interpstring(2:end),xxx);
        w = waitbar(0,sprintf('Measuring the voxel values (method %d) using %s interpolation for mouse %02d...',sampletype,interpstring(2:end),xxx));
        for i = 1:length(voxellist)
            waitbar(i/length(voxellist));
            pointindex(voxellist(i,1),voxellist(i,2)) = i;
            for q = -(step/2):(step/2)
                p = q+step/2+1;
                voxelvalue = (voxellist(i,:)+((q/f)*fvcorr(i,:)));
                x = voxelvalue(1);
                y = voxelvalue(2);
                % avoid going out of bounds or sampling excluded tissue
                % classes
                if (1 < x) && (x < size(MASK,1)-1) && (1 < y) && (y < size(MASK,2)-1)
                    if SEGS{xxx}(ceil(x),ceil(y)) > 0
                        for b = 1:numconts
                            sample_mat(i,p,b)  = interpn(data_proc{xxx}(:,:,b),x,y,interpstring);
                        end
                        traces(ceil(x),ceil(y)) = traces(ceil(x),ceil(y)) + 1;
                        lines(ceil(x),ceil(y)) = p;
                    else
                        break
                    end
                else
                    break
                end
            end
        end
    end
    
    close(w)
    
    % Save out some images to check profile positions and directionality
    TRA = make_nii(traces, pdims);
    save_nii(TRA,'Profiles.nii.gz');
    POI = make_nii(pointindex, pdims);
    save_nii(POI,'EdgeIndex.nii.gz');
    LIN = make_nii(lines, pdims);
    save_nii(LIN,'Positions.nii.gz');
    
    % Main profile properties
  
    sample_matfloat = double(sample_mat);
    fprintf('Calculating the properties of the normal profiles for mouse %02d...\n', xxx);
    slopec  = NaN(length(voxellist),numconts);
    tload   = NaN(length(voxellist),numconts);
    tloadc  = NaN(length(voxellist),numconts);
    fg      = NaN(length(voxellist),1);
    
    hg = waitbar(0,sprintf('Calculating the properties of the normal profiles for mouse %02d...',xxx));
    for j = 1:length(sample_matfloat)
        waitbar(j/length(sample_matfloat));
        for c = 1:numconts
            fg(j) = nnz(sample_mat(j,:,c));
            if fg(j) < 3*f
                tload(j,c) = 0;
            else
                tload(j,c) = sum(sample_matfloat(j,:,c));
                slopecx = (norm(fvcorr(j,:)):norm(fvcorr(j,:)):(step+1)*(norm(fvcorr(j,:))));
                slopecx = slopecx';
                slopecy = (nonzeros(sample_matfloat(j,:,c)));
                slopec(j,c) = slopecx(1:fg(j))\slopecy;
            end
        end
    end
    close(hg)
    clear j
    
    for hh = 1:numconts
        tloadc(:,hh) = tload(:,hh)./fg;
    end
    
    % stats
    meanADCprof = mean(sample_mat(:,:,adcnum),1);
    sdADCprof = std(sample_mat(:,:,adcnum),1);
    semADCprof = sdADCprof/sqrt(length(voxellist));
    meanHLAprof = mean(sample_mat(:,:,hlanum),1);
    sdHLAprof = std(sample_mat(:,:,hlanum),1);
    semHLAprof = sdHLAprof/sqrt(length(voxellist));
    meanHnEprof = mean(sample_mat(:,:,hnenum),1);
    sdHnEprof = std(sample_mat(:,:,hnenum),1);
    semHnEprof = sdHnEprof/sqrt(length(voxellist));
    
    if plotshow == 1 
     
        meanADCprof_ALL{xxx}=meanADCprof;
        meanHLAprof_ALL{xxx}=meanHLAprof;
        meanHnEprof_ALL{xxx}=meanHnEprof;
        semADCprof_ALL{xxx}=semADCprof;
        semHLAprof_ALL{xxx}=semHLAprof;
        semHnEprof_ALL{xxx}=semHnEprof;   
       
    end
    
    ADCprofiles = reshape(sample_mat(:,:,adcnum),size(sample_mat(:,:,adcnum),1)*size(sample_mat(:,:,adcnum),2),1);
    HLAprofiles = reshape(sample_mat(:,:,hlanum),size(sample_mat(:,:,hlanum),1)*size(sample_mat(:,:,hlanum),2),1);
    HnEprofiles = reshape(sample_mat(:,:,hnenum),size(sample_mat(:,:,hnenum),1)*size(sample_mat(:,:,hnenum),2),1);

    % MOAT scatterplot based on the widening tumour 'moat'
    ADC.img = data_proc{xxx}(:,:,adcnum);
    HLA.img = data_proc{xxx}(:,:,hlanum);
    HnE.img = data_proc{xxx}(:,:,hnenum);
    
    if scatshow == 1
        
        for ii = 1:moatsize
     
            ADCdata{xxx}=ADC.img(BW1wide{xxx}(:,:,ii)==1);
            HLAdata{xxx}=HLA.img(BW1wide{xxx}(:,:,ii)==1);

            A = HLAdata{xxx};
            B= ADCdata{xxx};
            [A,idx]=sort(A,'ascend');
            B=B(idx);
            
            if (xxx<11)
                A(A>0.95)=NaN
                A(A<0.05)=NaN
                B(B>0.85e-3)=NaN  % normal tissue 0.85e-3
                B(B<0.55e-3)=NaN  % 0.5 e-4
            end
            
            idx2 = isfinite(A) & isfinite(B);
            GbA{xxx}=A
            GbB{xxx}=B
            Gbidx{xxx}=idx2
            
            % standard ADC against SIH
            [f0, gof] = linFit(A(idx2), B(idx2))
            [RHO,PVAL] = corr(A(idx2), B(idx2),'Type','Spearman');
            figure;
            plot( f0,'--r', A(idx2), B(idx2),'o');
            axis([0.05 0.95 0.00055 0.00085])
            title(sprintf('Slope %.2d bsln %.2d R2 %.2d Rho %.2d pval %.2d',f0.p1, f0.p2, gof.rsquare,RHO,PVAL))
            saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHLAfitmouse%02dmoat%d',xxx,ii), '.tif'));
            saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHLAfitmouse%02dmoat%d',xxx,ii), '.pdf'));
            
            Invcoef{xxx}=f0.p1;
            Densnorm{xxx}=f0.p2;
            Fitrsquare{xxx}=gof.rsquare;
            GlobInvcoef{xxx}=f0.p1;
            GlobDensnorm{xxx}=f0.p2;
            GlobFitrsquare{xxx}=gof.rsquare;
           
            
            % standard SIH against ADC
            [f0, gof] = linFit(B(idx2), A(idx2))
            [RHO,PVAL] = corr(B(idx2), A(idx2),'Type','Spearman');
            figure;
            plot( f0,'--r', B(idx2), A(idx2),'o');
            axis([0.00055 0.00085 0.05 0.95 ])
            title(sprintf('Inv Slope %.2d bsln %.2d R2%.2d Rho %.2d pval %.2d' ,f0.p1, f0.p2, gof.rsquare,RHO,PVAL))
            saveas(gcf,strcat(outputfolder, sprintf('\\InvADCvHLAfitmouse%02dmoat%d',xxx,ii), '.tif'));
            saveas(gcf,strcat(outputfolder, sprintf('\\InvADCvHLAfitmouse%02dmoat%d',xxx,ii), '.pdf'));
       
            Invcoef2{xxx}=f0.p1;
            Densnorm2{xxx}=f0.p2;
            Fitrsquare2{xxx}=gof.rsquare;
            GlobInvcoef2{xxx}=f0.p1;
            GlobDensnorm2{xxx}=f0.p2;
            GlobFitrsquare2{xxx}=gof.rsquare;
            
            close all;

        end
        
    end
    
    moatADC = ADC.img(BW1wide{xxx}(:,:,moatsize)==1);
    moatHLA = HLA.img(BW1wide{xxx}(:,:,moatsize)==1);
    moatHnE = HnE.img(BW1wide{xxx}(:,:,moatsize)==1);
    
    % ztransform as ADC values are tiny and this may affect weighting of other
    % mathematical assessments such as clustering due to weightings...
    moatADCz = zscore(moatADC);
    moatHLAz = zscore(moatHLA);
    moatHnEz = zscore(moatHnE);
    
    % Make the whole cohort dataset
    if tu == 1
        
        indexvectorprof = ones(length(sample_matfloat),1);
        indexvectormoat = ones(length(moatADC),1);
        % THIS IS THE WHOLE DATASET FOR ALL MICE
        % SAMPLED FROM PROFILES
        fullprofiledata = sample_matfloat;
        % SAMPLED FROM THE BOUNDARY ROI
        fullADCdata = moatADC;
        fullHLAdata = moatHLA;
        fullHnEdata = moatHnE;
        % SAMPLED FROM THE BOUNDARY ROI AND Z TRANSFORMED
        fullADCdataz = moatADCz;
        fullHLAdataz = moatHLAz;
        fullHnEdataz = moatHnEz;
        
    else
        
        indexvectorprof = [indexvectorprof;xxx*ones(length(sample_matfloat),1)];
        indexvectormoat = [indexvectormoat;xxx*ones(length(moatADC),1)];
        % THIS IS THE WHOLE DATASET FOR ALL MICE
        % SAMPLED FROM PROFILES
        fullprofiledata = cat(1,fullprofiledata,sample_matfloat);
        % SAMPLED FROM THE BOUNDARY ROI
        fullADCdata = [fullADCdata;moatADC];
        fullHLAdata = [fullHLAdata;moatHLA];
        fullHnEdata = [fullHnEdata;moatHnE];
        % SAMPLED FROM THE BOUNDARY ROI AND Z TRANSFORMED
        fullADCdataz = [fullADCdataz;moatADCz];
        fullHLAdataz = [fullHLAdataz;moatHLAz];
        fullHnEdataz = [fullHnEdataz;moatHnEz];
        
    end
    
    tu = tu + 1;
    %CLUSTERING
    moatclust = [moatADC,moatHLA];
    moatclustz = [moatADCz,moatHLAz];
    
    clear C U index1 index2 index3
    [C,U] = fcm(moatclustz,3);
    maxU = max(U);
    index1 = find(U(1,:) == maxU);
    index2 = find(U(2,:) == maxU);
    index3 = find(U(3,:) == maxU);
    
    opts = statset('Display','final');
    [idx,C] = kmeans(moatclustz,3,'Distance','cityblock','Replicates',5,'Options',opts);

    idxF = zeros(length(moatclustz),1);
    idxF(index1) = 1;
    idxF(index2) = 2;
    idxF(index3) = 3;
    
    moatclustimzK = zeros(size(VOI));
    moatclustimzK(BW1wide{xxx}(:,:,moatsize)==1) = idx;
    MCZIK = make_nii(moatclustimzK, pdims);
    save_nii(MCZIK,'ClustClass_3KMCz_ADC_HLA.nii.gz');
    
    moatclustimzF = zeros(size(VOI));
    moatclustimzF(BW1wide{xxx}(:,:,moatsize)==1) = idxF;
    MCZIF = make_nii(moatclustimzF, pdims);
    save_nii(MCZIF,'ClustClass_3FCMz_ADC_HLA.nii.gz');
     
    % Correlation stats on the profiles between ADC and tissue measures
    rholist = NaN(1,length(voxellist));
    pvalist = NaN(1,length(voxellist));
    
    for uu = 1:length(voxellist)
        [RHO,PVAL] = corr([squeeze(sample_mat(uu,:,2));squeeze(sample_mat(uu,:,3))]', 'type', 'Spearman'); %QUESTION
        rholist(uu) = RHO(1,2);
        pvalist(uu) = PVAL(1,2);
        % Map for 1-slope p value UNCORRECTED!
        pvalscompnc(voxellist(uu,1),voxellist(uu,2)) = 1 - pvalist(uu);
        % Correction (BONFERRONI )
        pvalscompc(voxellist(uu,1),voxellist(uu,2)) = 1 - (pvalist(uu).*length(voxellist));
        % Actual correlation coefficient (moment, rho, tau, depending)
        corrcoef(voxellist(uu,1),voxellist(uu,2)) = rholist(uu);
    end
    
    tc_ADC = zeros(size(VOI));
    tc_HLA = zeros(size(VOI));
    tc_HnE = zeros(size(VOI));
    for uu = 1:length(voxellist)
        tc_ADC(voxellist(uu,1),voxellist(uu,2)) = slopec(uu,adcnum);
        tc_HLA(voxellist(uu,1),voxellist(uu,2)) = slopec(uu,hlanum);
        tc_HnE(voxellist(uu,1),voxellist(uu,2)) = slopec(uu,hnenum);
    end
    
    % save out correlation coefficient map for inspection
    COCO = make_nii(corrcoef, pdims);
    save_nii(COCO,'Correlation_coef.nii.gz');
    
    % save out p-value map for visual inspection
    % keep alpha at 0.05 for visualisation by effectively multiplying the pvalue by the number of nulls
    PVNC = make_nii(pvalscompnc, pdims);
    save_nii(PVNC,'1minusPvals_nc.nii.gz');
    
    PVC = make_nii(pvalscompc, pdims);
    save_nii(PVNC,'1minusPvals_bc.nii.gz');
    
    pcsig = 100*(length(find((pvalist.*length(voxellist))<0.05))/length(voxellist));
    
    % save out slope maps for inspection
    TC_ADC = make_nii(tc_ADC, pdims);
    save_nii(TC_ADC,'TC_ADC.nii.gz');
    TC_HLA = make_nii(tc_HLA, pdims);
    save_nii(TC_HLA,'TC_HLA.nii.gz');
    TC_HnE = make_nii(tc_HnE, pdims);
    save_nii(TC_HnE,'TC_HnE.nii.gz');
    
    if boxpshow == 1
         
    end
    
    voxellist_ALL {xxx}=voxellist;
    ft_ALL{xxx}=ft;
    cd ..
    
end


%% Checks and outputs
%% PLOTS

%% All ADC vs HLA allscatter
close all;
figure;
for xxx = mouselist
    
    if (xxx~=9)||(xxx~=14)
        scatter(HLA.img(BW1wide{xxx}(:,:,ii)==1),ADC.img(BW1wide{xxx}(:,:,ii)==1),'k o');
        hold on
    end
    
end
% axis ( [0.1 0.9 0.4e-3 0.9e-3])
title(sprintf('Scatter ADC vs HLA all'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHLAscater'), '.tif'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHLAscater'), '.pdf'));
%% All ADC vs HnE scatter
close all;
figure;
for xxx = mouselist
    
    if (xxx~=3)||(xxx~=7)||(xxx~=8)||(xxx~=9)||(xxx~=14)
        scatter(HnE.img(BW1wide{xxx}(:,:,ii)==1),ADC.img(BW1wide{xxx}(:,:,ii)==1));
        hold on
    end
    
end
% axis ( [0.2 1 0.4e-3 0.9e-3])
title(sprintf('Scatter ADC vs H&E all'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHnEscater'), '.tif'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHnEscater'), '.pdf'));


%% Vector data: ADC profile, HLA profile, ADC vs HLA profile 
close all;
for xxx = mouselist
    
    a=0.6;
    b=0;
    c=0;
    color=[a 0 0];
    
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    col='-k'

    
    switch xxx
        case 11; col='-r'
        case 12; col='-k'
        case 13; col='-b'
    end
    % ADC vs Histology
    figure('Name',sprintf('ADCvsHLA - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [310 10 1400 500],'color','w')
    subtightplot(1,3,1,0.05);
    plot(meanADCprof_ALL{xxx},col);
    hold on
    errorbar(meanADCprof_ALL{xxx},semADCprof_ALL{xxx},col);
%     ylim( [0.55e-3 0.75e-3])
xlim([0 14 ]);
    title(sprintf('Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,2,0.05);
    plot(meanHLAprof_ALL{xxx},col);
    hold on
    errorbar(meanHLAprof_ALL{xxx},semHLAprof_ALL{xxx},col);
% ylim( [0.1 0.9])
xlim([0 14 ]);
    title(sprintf('Mean HLA Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,3,0.05);
    plot(meanHLAprof_ALL{xxx},meanADCprof_ALL{xxx},col);
    hold on
     errorbarxy(meanHLAprof_ALL{xxx}, meanADCprof_ALL{xxx}, semHLAprof_ALL{xxx}, semADCprof_ALL{xxx},col);
    
% axis( [0.1 0.9 0.55e-3 0.75e-3])
    title(sprintf('Mean HLA vs Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
    saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto%d',xxx), '.pdf'));
    
    %
    
end
%% %Plot all together ADC/HLA
close all;

figure('Name',sprintf('ADCvsHLA_model'),'NumberTitle','off');
for xxx = mouselist
    
    a=0.6;
    b=0;
    c=0;
    color=[a 0 0];
    switch xxx
        case 11; col='-r'
        case 12; col='-k'
        case 13; col='-b'
    end
    
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,1)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    % ADC vs Histology
    
    set(gcf, 'Position', [310 10 1400 500],'color','w')
    subtightplot(1,3,1,0.05);
    
    plot(meanADCprof_ALL{xxx},col);
    hold on
    errorbar(meanADCprof_ALL{xxx},semADCprof_ALL{xxx},col);
%     ylim( [0.3e-3 0.7e-3])
    title(sprintf('Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,2,0.05);
    plot(meanHLAprof_ALL{xxx},col);
    hold on
    errorbar(meanHLAprof_ALL{xxx},semHLAprof_ALL{xxx},col);
    %     ylim( [0.1 0.9])
    title(sprintf('Mean HLA Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,3,0.05);
    plot(meanHLAprof_ALL{xxx},meanADCprof_ALL{xxx},col);
    hold on
    errorbarxy(meanHLAprof_ALL{xxx}, meanADCprof_ALL{xxx}, semHLAprof_ALL{xxx}, semADCprof_ALL{xxx},col);
    
end
%     axis( [0.1 0.9 0.4e-3 1e-3])
title(sprintf('Mean HLA vs Mean ADC Profile with SEM Error Bars'))
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto_model'), '.tif'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto_model'), '.pdf'));

%% %Plot all together ADC/HnE
close all;

figure('Name',sprintf('ADCvsHnE_model'),'NumberTitle','off');
for xxx = mouselist
    
    a=0.6;
    b=0;
    c=0;
    color=[a 0 0];
    switch xxx
        case 11; col='-r'
        case 12; col='-k'
        case 13; col='-b'
    end
    
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,1)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    % ADC vs Histology
    
    set(gcf, 'Position', [310 10 1400 500],'color','w')
    subtightplot(1,3,1,0.05);
    
    plot(meanADCprof_ALL{xxx},col);
    hold on
    errorbar(meanADCprof_ALL{xxx},semADCprof_ALL{xxx},col);
%     ylim( [0.3e-3 0.7e-3])
    title(sprintf('Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,2,0.05);
    plot(meanHnEprof_ALL{xxx},col);
    hold on
    errorbar(meanHnEprof_ALL{xxx},semHnEprof_ALL{xxx},col);
    %     ylim( [0.1 0.9])
    title(sprintf('Mean HnE Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,3,0.05);
    plot(meanHnEprof_ALL{xxx},meanADCprof_ALL{xxx},col);
    hold on
    errorbarxy(meanHnEprof_ALL{xxx}, meanADCprof_ALL{xxx}, semHnEprof_ALL{xxx}, semADCprof_ALL{xxx},col);
    
end
%     axis( [0.1 0.9 0.4e-3 1e-3])
title(sprintf('Mean HnE vs Mean ADC Profile with SEM Error Bars'))
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto_model'), '.tif'));
saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHisto_model'), '.pdf'));

%% Plot HLA v HnE
close all;

for xxx = mouselist
    % Histology
    figure('Name',sprintf('Histology - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [0 0 1400 500],'color','w')
    subtightplot(1,3,1,0.05);
    plot(meanHnEprof_ALL{xxx},'-k');
    hold on
    errorbar(meanHnEprof_ALL{xxx},semHnEprof_ALL{xxx},'-k');
    title(sprintf('Mean H&E Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,2,0.05);
    plot(meanHLAprof_ALL{xxx},'-k');
    hold on
    errorbar(meanHLAprof_ALL{xxx},semHLAprof_ALL{xxx},'-k');
    title(sprintf('Mean HLA Profile with SEM Error Bars Mouse %02d',xxx))
    
    subtightplot(1,3,3,0.05);
    plot(meanHLAprof_ALL{xxx},meanHnEprof_ALL{xxx},'-k');
    hold on
    errorbarxy(meanHLAprof_ALL{xxx}, meanHnEprof_ALL{xxx}, semHLAprof_ALL{xxx}, semHnEprof_ALL{xxx},'--r');
    title(sprintf('Mean HLA vs Mean HnE Profile with SEM Error Bars Mouse %02d',xxx))
    saveas(gcf,strcat(outputfolder, sprintf('\\Histo%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\Histo%d',xxx), '.pdf'));
    
end

%% Plot ADC wo/with Vector overlays

for xxx = mouselist
    % Vectors
    figure('Name',sprintf('T2W Vector overlays - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [0 0 800 400],'color','w')
    subtightplot(1,2,1,0.01);
    imagesc( T2Wim{xxx},[0.5 1] ); colormap(gray);
    subtightplot(1,2,2,0.01);
    imagesc( T2Wim{xxx} ); colormap(gray);
    hold on ; quiver(voxellist_ALL{xxx}(:,2),voxellist_ALL{xxx}(:,1),ft_ALL{xxx}(:,2),ft_ALL{xxx}(:,1),2, '-g');
    saveas(gcf,strcat(outputfolder, sprintf('\\VectT2W%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\VectT2W%d',xxx), '.pdf'));
    
    
    figure('Name',sprintf('ADC Vector overlays - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [0 0 800 400],'color','w')
    subtightplot(1,2,1,0.01);
    imagesc( ADCim{xxx},[0.2 0.6] ); colormap(gray);
    subtightplot(1,2,2,0.01);
    imagesc( ADCim{xxx} ); colormap(gray);
    hold on ; quiver(voxellist_ALL{xxx}(:,2),voxellist_ALL{xxx}(:,1),ft_ALL{xxx}(:,2),ft_ALL{xxx}(:,1),2, '-g');
    %     hold on ; quiver(voxellist_ALL{xxx}(:,2),voxellist_ALL{xxx}(:,1),fvcorr(:,2),fvcorr(:,1),2, '-g');
    
    saveas(gcf,strcat(outputfolder, sprintf('\\VectADC%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\VectADC%d',xxx), '.pdf'));
    
end
%% Plot MOAT masks overlays

for xxx = mouselist
    
    % MOATS
    moat2ovelay{xxx} = imoverlay(ADCim{xxx}, mat2gray(squeeze(BW1wide{xxx}(:,:,2))), color);
    moat4ovelay{xxx} = imoverlay(ADCim{xxx}, mat2gray(squeeze(BW1wide{xxx}(:,:,4))), color);
    moat6ovelay{xxx} = imoverlay(ADCim{xxx}, mat2gray(squeeze(BW1wide{xxx}(:,:,6))), color);
    moat8ovelay{xxx} = imoverlay(ADCim{xxx}, mat2gray(squeeze(BW1wide{xxx}(:,:,8))), color);
    
    figure('Name',sprintf('MOAT overlays 2 4 6 8 - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [0 0 800 800],'color','w')
    subtightplot(2,2,1,0.01);
    imagesc(moat2ovelay{xxx}); colormap(gray);
    subtightplot(2,2,2,0.01);
    imagesc(moat4ovelay{xxx});colormap(gray);
    subtightplot(2,2,3,0.01);
    imagesc(moat6ovelay{xxx});colormap(gray);
    subtightplot(2,2,4,0.01);
    imagesc(moat8ovelay{xxx});colormap(gray);
    saveas(gcf,strcat(outputfolder, sprintf('\\Moat%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\Moat%d',xxx), '.pdf'));
    
end


%% Supp Info Summary figures

for xxx = mouselist

    fontSize=14
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    col='-k'
    
    moat2ovelay{xxx} = imoverlay(ADCim{xxx}, mat2gray(squeeze(BW1wide{xxx}(:,:,2))), color);
    
    % ADC vs Histology
    figure('Name',sprintf('ADCvsHLA - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [10 10 1700 400],'color','w')
    
    marg_h=0.1
    marg_w=0.02
    subtightplot(1,5,1,0.03,marg_h,marg_w);
    imagesc(moat2ovelay{xxx}); colormap(gray);
    title(sprintf('ADC - ROI M%02d',xxx),'FontSize', fontSize)
    
    subtightplot(1,5,2,0.03,marg_h,marg_w);
    imagesc( ADCim{xxx} ); colormap(gray);
    hold on ; quiver(voxellist_ALL{xxx}(:,2),voxellist_ALL{xxx}(:,1),ft_ALL{xxx}(:,2),ft_ALL{xxx}(:,1),2, '-g');
    title(sprintf('ADC - Vect Profiles M%02d',xxx),'FontSize', fontSize)
  
    subtightplot(1,5,3,0.05,marg_h,marg_w);
    plot(meanADCprof_ALL{xxx},col);
    hold on
    errorbar(meanADCprof_ALL{xxx},semADCprof_ALL{xxx},col);
    %     ylim( [0.55e-3 0.75e-3])
    xlim([0 14 ]);
    title(sprintf('Mean ADC Prof. M%02d',xxx),'FontSize', fontSize)
    xlabel('profile point') 
    ylabel('ADC / mm s-1') 
 
    subtightplot(1,5,4,0.05,marg_h,marg_w);
    plot(meanHLAprof_ALL{xxx},col);
    hold on
    errorbar(meanHLAprof_ALL{xxx},semHLAprof_ALL{xxx},col);
    xlim([0 14 ]);
    title(sprintf('Mean HLA Prof. M%02d',xxx),'FontSize', fontSize)
    xlabel('profile point') 
    ylabel('HLA / a.u.') 
    
    subtightplot(1,5,5,0.05,marg_h,marg_w);
    plot(meanHLAprof_ALL{xxx},meanADCprof_ALL{xxx},col);
    hold on
    errorbarxy(meanHLAprof_ALL{xxx}, meanADCprof_ALL{xxx}, semHLAprof_ALL{xxx}, semADCprof_ALL{xxx},col);
    title(sprintf('HLA / ADC Prof. M%02d',xxx),'FontSize', fontSize)
    xlabel('HLA / a.u.') 
    ylabel('ADC / mm s-1') 
    
    saveas(gcf,strcat(outputfolder, sprintf('\\SuppFigAll%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\SuppFigall%d',xxx), '.pdf'));
      
       
end



%% ROIS 

close all;

ADCclims{1}(1,1)=0.3;
ADCclims{1}(1,2)=0.8;
T2clims{1}(1,1)=0.5;
T2clims{1}(1,2)=0.9;

ADCclims{3}(1,1)=0.4;
ADCclims{3}(1,2)=0.8;
T2clims{3}(1,1)=0.5;
T2clims{3}(1,2)=0.9;

ADCclims{4}(1,1)=0.3;
ADCclims{4}(1,2)=0.9;
T2clims{4}(1,1)=0.6;
T2clims{4}(1,2)=0.9;

ADCclims{5}(1,1)=0.5;
ADCclims{5}(1,2)=0.8;
T2clims{5}(1,1)=0.5;
T2clims{5}(1,2)=0.9;

ADCclims{6}(1,1)=0.5;
ADCclims{6}(1,2)=0.9;
T2clims{6}(1,1)=0.6;
T2clims{6}(1,2)=0.9;

ADCclims{7}(1,1)=0.4;
ADCclims{7}(1,2)=0.8;
T2clims{7}(1,1)=0.5;
T2clims{7}(1,2)=1;

ADCclims{8}(1,1)=0.5;
ADCclims{8}(1,2)=0.9;
T2clims{8}(1,1)=0.5;
T2clims{8}(1,2)=1;

ADCclims{9}(1,1)=0.4;
ADCclims{9}(1,2)=0.8;
T2clims{9}(1,1)=0.5;
T2clims{9}(1,2)=0.9;

ADCclims{10}(1,1)=0.4;
ADCclims{10}(1,2)=0.7;
T2clims{10}(1,1)=0.6;
T2clims{10}(1,2)=0.9;

ADCclims{12}(1,1)=0.4;
ADCclims{12}(1,2)=0.8;
T2clims{12}(1,1)=0.5;
T2clims{12}(1,2)=0.9;

ADCclims{11}(1,1)=0.4;
ADCclims{11}(1,2)=0.8;
T2clims{11}(1,1)=0.5;
T2clims{11}(1,2)=0.9;

ADCclims{13}(1,1)=0.4;
ADCclims{13}(1,2)=0.8;
T2clims{13}(1,1)=0.5;
T2clims{13}(1,2)=0.9;

ADCclims{14}(1,1)=0.4;
ADCclims{14}(1,2)=0.8;
T2clims{14}(1,1)=0.5;
T2clims{14}(1,2)=0.9;

for xxx = mouselist
    
    %ROIs
    a=0.6;
    b=0;
    c=0;
    color=[a 0 0];
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    HLAim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,3)));
    if (xxx==9)
        HLAim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,4)));
    end
    
    HnEm{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,4)));
    
    Roiovelay{xxx} = imoverlay(T2Wim{xxx}, MASKS{xxx}, color);
    cont1ovelay{xxx} = imoverlay(T2Wim{xxx}, BW1{xxx}, color);
    cont2ovelay{xxx} = imoverlay(T2Wim{xxx}, BW2{xxx}, color);
    
    
    ADCRoiovelay{xxx} = imoverlay(ADCim{xxx}, MASKS{xxx}, color);
    T2WRoiovelay{xxx} = imoverlay(T2Wim{xxx}, T2WMASKS{xxx}, color);
    HLARoiovelay{xxx} = imoverlay(HLAim{xxx}, HLAMASKS{xxx}, color);
    
    figure('Name',sprintf('ROI & contour overlays - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [0 0 800 800],'color','w')
    subtightplot(2,2,1,0.01);
    imagesc(squeeze(T2Wim{xxx}),[0.5 1]); colormap(gray);
    subtightplot(2,2,2,0.01);
    imagesc(Roiovelay{xxx});colormap(gray);
    subtightplot(2,2,3,0.01);
    imagesc(cont1ovelay{xxx});colormap(gray);
    subtightplot(2,2,4,0.01);
    imagesc(cont2ovelay{xxx});colormap(gray);
    saveas(gcf,strcat(outputfolder, sprintf('\\ROIs%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\ROIs%d',xxx), '.pdf'));
    
    
    figure('Name',sprintf('ROI & contour overlays - mouse %d', xxx),'NumberTitle','off');
    set(gcf, 'Position', [110 100 800 800],'color','w')
    
    
    subtightplot(3,2,1,0.01);
    imagesc(ADCim{xxx},[ADCclims{xxx}(1,1) ADCclims{xxx}(1,2)]); colormap(gray);
    subtightplot(3,2,2,0.01);
    imagesc(ADCRoiovelay{xxx});colormap(gray);
    subtightplot(3,2,3,0.01);
    imagesc(squeeze(T2Wim{xxx}),[T2clims{xxx}(1,1) T2clims{xxx}(1,2)]); colormap(gray);
    subtightplot(3,2,4,0.01);
    imagesc(T2WRoiovelay{xxx});colormap(gray);
    subtightplot(3,2,5,0.01);
    imagesc(HLAim{xxx});colormap(gray);
    subtightplot(3,2,6,0.01);
    imagesc(HLARoiovelay{xxx});colormap(gray);
    saveas(gcf,strcat(outputfolder, sprintf('\\ROIsv2%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\ROIsv2%d',xxx), '.pdf'));
    
end

%% Vector fits

for xxx = mouselist
    xxx
    a=0.6;b=0; c=0;
    color=[a 0 0];
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    meanHLApnorm=meanHLAprof_ALL{xxx};
    meanADCpnorm=meanADCprof_ALL{xxx};
    
    [f1HLA, HLAgof] = linFit(meanHLApnorm,meanADCpnorm)
    Invcoef{xxx}=f1HLA.p1;
    Densnorm{xxx}=f1HLA.p2;
    Fitrsquare{xxx}=HLAgof.rsquare;
    
    figure;
    plot( f1HLA,'--r', meanHLApnorm,meanADCpnorm,'o');
    %     axis([0.1 0.9 0.4e-3 0.85e-3])
    title(sprintf('ADC vs HLA lin fit (a=%03d,  b= %03d R2= %03d) Mouse %02d',Invcoef{xxx},Densnorm{xxx},  Fitrsquare{xxx},xxx))
    saveas(gcf,strcat(outputfolder, sprintf('\\ADCHLvectAfit%d',xxx), '.tif'));
    saveas(gcf,strcat(outputfolder, sprintf('\\ADCHLvectAfit%d',xxx), '.pdf'));
    
end


%% ADC v HnE plots 
close all;
for xxx = mouselist
    
    a=0.6;
    b=0;
    c=0;
    color=[a 0 0];
    T2Wim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    ADCim{xxx}=mat2gray(squeeze(data_proc{xxx}(:,:,2)));
    
    
    if (xxx~=3)||(xxx~=7)||(xxx~=8)||(xxx~=11)||(xxx~=12)||(xxx~=13)
        % ADC vs Histology
        figure('Name',sprintf('ADCvsHnE - mouse %d', xxx),'NumberTitle','off');
        set(gcf, 'Position', [310 10 1400 500],'color','w')
        subtightplot(1,3,1,0.05);
        plot(meanADCprof_ALL{xxx},'-k');
        hold on
        errorbar(meanADCprof_ALL{xxx},semADCprof_ALL{xxx},'-k');
        title(sprintf('Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
        
        subtightplot(1,3,2,0.05);
        plot(meanHnEprof_ALL{xxx},'-k');
        hold on
        errorbar(meanHnEprof_ALL{xxx},semHnEprof_ALL{xxx},'-k');
        title(sprintf('Mean HnE Profile with SEM Error Bars Mouse %02d',xxx))
        
        subtightplot(1,3,3,0.05);
        plot(meanHnEprof_ALL{xxx},meanADCprof_ALL{xxx},'-k');
        hold on
        errorbarxy(meanHnEprof_ALL{xxx}, meanADCprof_ALL{xxx}, semHnEprof_ALL{xxx}, semADCprof_ALL{xxx},'--r');
        title(sprintf('Mean HnE vs Mean ADC Profile with SEM Error Bars Mouse %02d',xxx))
        saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHnE%d',xxx), '.tif'));
        saveas(gcf,strcat(outputfolder, sprintf('\\ADCvHnE%d',xxx), '.pdf'));
        
    end
    
    
end
%% Create Results matrix
for xxx = mouselist
    
    ADCpixvol{xxx}=sum(sum(MASKS{xxx}))
    T2Wpixvol{xxx}=sum(sum(T2WMASKS{xxx}))
    HLApixvol{xxx}=sum(sum(HLAMASKS{xxx}))
    HLAoverhead{xxx}=(HLApixvol{xxx}-ADCpixvol{xxx})/ADCpixvol{xxx}
    T2Woverhead{xxx}=(T2Wpixvol{xxx}-ADCpixvol{xxx})/ADCpixvol{xxx}
    HLARatio{xxx}=HLApixvol{xxx}/ADCpixvol{xxx}
    T2WRatio{xxx}=T2Wpixvol{xxx}/ADCpixvol{xxx}
    Invcoef{xxx}=Invcoef{xxx}
    Densnorm{xxx}=Densnorm{xxx}
    Fitrsquare{xxx}=Fitrsquare{xxx}
    
end


i=1
for xxx = mouselist
    
    tResults(i,1)=ADCpixvol{xxx}
    tResults(i,2)=HLApixvol{xxx}
    tResults(i,3)=HLAoverhead{xxx}
    tResults(i,4)=HLARatio{xxx}
    tResults(i,5)=Invcoef{xxx}
    tResults(i,6)=Densnorm{xxx}
    tResults(i,7)=Fitrsquare{xxx}
    tResults(i,8)=T2Wpixvol{xxx}
    tResults(i,9)=T2Woverhead{xxx}
    tResults(i,10)=T2WRatio{xxx}
    tResults(i,11)=xxx
    
    i=i+1
end


%% Save Workspace

%  Save Workspace
save(strcat(outputfolder, sprintf('\\workspace.mat'))); % Type in name of file.

