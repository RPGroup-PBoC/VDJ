%[X,Y,fps,success,fpw,clA]=loadrawV3(path, filenames, varargin)
%
%Steph's best version of code that loads raw TPM data (the POS files that
%are the output of realtime).  See open_pxlV3 for code that reads the pxl
%movie files.
%
%Inputs:
%path: full path of the folder which contains the folder of data to load
%filename: name of the folder to load (and the first part of all the files
%in the folder)
%Optional input: "Nolac1" or "Lac2" or something that specifies for the
%user which folder is attempting to be loaded (in case there's an error in
%the loading process when this is called from masterscriptV3)
%
%Outputs: X and Y are bead-by-frames matrices containing X and Y data for
%each bead respectively.  fps is the frame rate of the camera for this data
%set.  Success is 1 if the data were loaded successfully, and 0 if not.
%Fpw is the number of frames per file (used by masterscriptV3 for loading
%pxl files to watch movies during analysis).  clA is also used by
%masterscriptV3 and is the integer class of the pxl files (with our current
%version of acquisition code, this is uint16).
%
%Stephanie Johnson 1/11

function [X,Y,fps,success,fpw,clA]=loadrawV3(path, filename, varargin)

defaultclA = 'uint16';

success = 1;

loadname = fullfile(path,filename,filename);

if nargin == 2
    name = 'this set';
else
    name = varargin{1};
end

%Load the ANAL file, which is the only file that keeps track of how many
%total files are in the folder (that got analyzed).  It also contains the
%number of beads, though INFO and each file also contains that.
if exist(strcat(loadname,'ANAL.mat')) %If stop wasn't run, the ANAL file won't be created, so need to check
    prelim = load(strcat(loadname, 'ANAL.mat'));
    %bds = size(prelim.rms_pos,1);%How many beads in this set %***
    numfiles = prelim.win_anal;
else
    usrinput = input(strcat('Cannot find ANAL.mat for ',name,'; continue? (y/n)'),'s');
    if strcmpi(usrinput,'y')
        numfiles = input(strcat('How many POS files in directory ',fullfile(path,filenames),'?' ));
        firstfile = load(strcat(loadname,'_1_POS'));
        bds(i)=firstfile.nroi;
    elseif strcmpi(usrinput,'n')
        success=0;
        X=0;
        Y=0;
        return
    end
end

infofile = load(strcat(loadname, 'INFO.mat')); %Unlike the ANAL file I think the INFO file is always created.
        %This file contains all the information about the imaging
        %conditions.
bds = size(infofile.roi,1);%***
fpw = infofile.fpw; %This is how many frames were stored in each file, need to load them in the code below
fps = 1/(infofile.acc_time); %Frame rate at which these data were taken (same as kin_time)
if isfield(infofile,'clA') %Older version of the acquisition code didn't save this
    clA = infofile.clA;
else
    %All our code uses uint16; uncomment this to allow user input to change
    %this default:
%     usrclA = input(strcat('Use clA = uint16 for',name,'? (y/n)'),'s');
%     if strcmpi(usrclA,'y')
%         clA = 'uint16';
%     else
%         clA = input('Use what clA?: ','s');
%     end
    clA = defaultclA;
end
    
tempx = zeros(numfiles, fpw, bds); %data is stored in each _POS.mat file as [frame (1 to wind),bd,x or y];
                                  %want to store it as [4-second chunk, frames in that chunk,
                                  %beads] and then reshape it below
tempy = zeros(numfiles, fpw, bds);
for j = 1:numfiles %Read in the raw x-y data from the pm_xy1 variable in each _POS.mat file
    data = load(strcat(loadname, '_' ,int2str(j), '_POS.mat')); 
    tempx(j, :, :)= data.pm_xy1(:, : , 1);
    tempy(j, :, :)= data.pm_xy1(:, : , 2);
    %strcat('Loaded file ',int2str(j),'/',int2str(numfiles(i)))   
    clear data
end

%Reshape into the desired beads-by-frames output matrix.  There's probably
%a way to avoid the for-loop but it wouldn't be a huge speed improvement
%(the number of beads is always relatively small) and reshape is
%confusing enough to use on 2-d matrices.

X = zeros(bds,numfiles*fpw);
Y = zeros(bds,numfiles*fpw);

for r=1:bds
   X(r,:) = reshape(transpose(tempx(:,:,r)),numfiles*fpw,1);
   Y(r,:) = reshape(transpose(tempy(:,:,r)),numfiles*fpw,1);
end

%When beads break, their x,y positions sometimes get set to NaN.  If you 
%fft something with a NaN, as in the drift correction and GaussFilt code,
%the whole thing becomes a NaN.  So set all NaNs to zeros:

X(isnan(X))=0;
Y(isnan(Y))=0;
    
