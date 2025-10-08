
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for preprossessing co-registered mouse data
%
% Copyright
% Dr Gerard Thompson, Dr Antoine Vallatos
% University of Edinburgh
% gerard.thompson@ed.ac.uk - antoine.vallatos@glasgow.ac.uk
%
% Version 1
% Tuesday 15th Oct 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
set(0,'defaultfigurecolor',[1 1 1])
%
% %%  Define input output folders
% Run and select the mouse folder (go in the data folder)
uiwait(helpdlg('Select input folder',...
    ' data analysis'));
inputfolder = uigetdir

model=1 %0 for data 1 for simulations

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

%populating the data matrix (dert)
saveoutslices = 0;
saveout = 1;
loadall = 1; %If this condition is true, it means that there is a requirement to load all data.
loaddert = 0; %If this condition is true, it means that the dert variable has not been loaded from any external source.
savedert = 1; %If this condition is true, it means that some data has been saved in the dert variable.
pdims = [0.11 0.11 1.5]; %depends on imaging but already harmonised on this dataset

% This is set to populate a multispectral sampling result matrix
numconts = 16;
plotshow = 1;
scatshow = 1;
clusshow = 1;
boxpshow = 1;

% Data provided here as a 3D multispectral stack of a single slice
% 1.T1w 2.T2w 3.T2map 4. T2whigh res1 5. T2whigh res2 6. T2whigh res3 7. FA 8. DWI 9. ADC 10. mbASLG 16%
% 11. T1Gd 12. mbASLG0% 13. HLA (M02 -> no HLA -  M06 HLA data not scaled) 14.	AD 15. RD 16. H&E. (M03, M07, M08 3,7,8 no H&E).

% For figures etc...
t1wnum = 1;
t2wnum = 2;
t2rnum = 3;
t2anum = 4;
t2bnum = 5;
t2cnum = 6;
fadnum = 7;
dwinum = 8;
adcnum = 9;
a16num = 10;
t1gnum = 11;
a00num = 12;
hlanum = 13;
addnum = 14;
rddnum = 15;
hnenum = 16;

labelvec = ['T1w','T2w','T2','T2wHR 1','T2wHR 2','T2wHR 3','FA','DWI','ADC','mbASLG16','T1wGd','mbASLG0','HLA','AD','RD','H&E'];
%% 

% Initialising dert as cell array with one element
if savedert == 1 && loaddert == 0 && loadall == 1
    dert = cell(1);
    
    %% Prepare normalised data mat files for processing
    
    for xxx = 1:10
        xxx
        load(sprintf('matlabM%02dwk12.mat',xxx));
        dert{xxx} = volume_image;
         if size(volume_image,3) < numconts
             dert{xxx}(:,:,size(volume_image,3)+1:16) = zeros(size(volume_image,1),size(volume_image,2));
         end
        if saveoutslices == 1
            hco = make_nii(dert{xxx});
            dirname = sprintf('M%02dwk12',xxx);
            mkdir(dirname);
            cd(dirname);
            for yyy = 1:numconts
                hcp = make_nii(hco.img(:,:,yyy));
                filename = sprintf('M%02dwk12_%02d.nii.gz',xxx,yyy);
                save_nii(hcp,filename);
            end
            cd ..
            filename2 = sprintf('M%02dwk12.nii.gz',xxx);
            save_nii(hco,filename2);
        end
    end
    
    
 
    for xxx = 11:13
        xxx
        load(sprintf('matlabM%02dwk12.mat',xxx));
        dert{xxx} = volume_image;
        if saveoutslices == 1
            hco = make_nii(dert{xxx});
            dirname = sprintf('M%02dwk12',xxx);
            mkdir(dirname);
            cd(dirname);
            for yyy = 1:numconts
                hcp = make_nii(hco.img(:,:,yyy));
                filename = sprintf('M%02dwk12_%02d.nii.gz',xxx,yyy);
                save_nii(hcp,filename);
            end
            cd ..
            filename2 = sprintf('M%02dwk12.nii.gz',xxx);
            save_nii(hco,filename2);
        end
    end
    
    save('DERT_raw.mat','dert');
    
%% Data corrections and adjustments
elseif savedert == 0 && loaddert == 1 && loadall == 0
    load('DERT_raw.mat');
else
    fprintf('Incorrect setting for data handling...')
end

lo_out = 0.5e-3
hi_out = 0.85e-3
range_out = hi_out - lo_out

for xxx=11:13
%Normalize to range [0, 1]:
dert{xxx}(:,:,2) = (dert{xxx}(:,:,2) - min(min(dert{xxx}(:,:,2)))) ./ (max(max(dert{xxx}(:,:,2)))- min(min(dert{xxx}(:,:,2))))
% %Normalize to range [lo_out, hi_out]:
lo_in = min(min(dert{xxx}(:,:,2))); %Minimum of each row
hi_in = max(max(dert{xxx}(:,:,2))); %Maximum of each row
range_in = hi_in - lo_in; %Range of each row
dert{xxx}(:,:,2) = ((dert{xxx}(:,:,2) - lo_in) ./ range_in) * range_out + lo_out
dert{xxx}(:,:,1)=dert{xxx}(:,:,2);

%Normalize to range [0, 1]:
dert{xxx}(:,:,3) = (dert{xxx}(:,:,3) - min(min(dert{xxx}(:,:,3)))) ./ (max(max(dert{xxx}(:,:,3)))- min(min(dert{xxx}(:,:,3))))
dert{xxx}(:,:,4)=dert{xxx}(:,:,3);
end





%Normalise m6 histology =(value-min)/(max-min)
dert{6}(:,:,13)=(dert{6}(:,:,13)-min(min(dert{6}(:,:,13))))./(max(max(dert{6}(:,:,13)))-min(min(dert{6}(:,:,13))));

%Normalise HnE (normalise=(value-min)/(max-min))
for i=1:10
    i
    dert{i}(:,:,16)=(dert{i}(:,:,16)-min(min(dert{i}(:,:,16))))./(max(max(dert{i}(:,:,16)))-min(min(dert{i}(:,:,16))));
end

save('DERT.mat','dert');


%% Dataset for analysis

for xxx = 1:10
data_proc{xxx}(:,:,1)=dert{xxx}(:,:,2);
data_proc{xxx}(:,:,2)=dert{xxx}(:,:,9);
data_proc{xxx}(:,:,3)=dert{xxx}(:,:,13);
data_proc{xxx}(:,:,4)=dert{xxx}(:,:,16);
end

 for xxx = 11:13
data_proc{xxx}(:,:,1)=dert{xxx}(:,:,1);
data_proc{xxx}(:,:,2)=dert{xxx}(:,:,2);
data_proc{xxx}(:,:,3)=dert{xxx}(:,:,3);
data_proc{xxx}(:,:,4)=dert{xxx}(:,:,4);
 end

 save('data_proc.mat','data_proc');
