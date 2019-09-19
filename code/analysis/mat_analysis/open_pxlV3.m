%function open_pxlV3(pathread, filename)
%
%Function to read a single .pxl file from data stored by realtime.  Inputs
%are the full path to the folder that contains the pxl files, and the name
%of the file to open (which is the last part of pathread followed 
%by _<#>.pxl).  This plays the entire .pxl file and is called by 
%masterscriptV3.  See PXLtomatrix to load the data in a pxl file into a
%matrix.
%
%Stephanie Johnson 1/11

function d = open_pxlV3(pathread, filename, varargin)

if nargin<3
    clA='uint16'; %Default class of data in the pxl files is uint16
else
    clA = varargin{1};
end

file=fullfile(pathread,filename);

fp = fopen(file,'r');%Open file read-only
fpw = fread(fp,1,clA); %Reads one element from the file fp, sets pointer at the second; first element must be ... ? number of frames per file I think
nroi = fread(fp,1,clA);%Second element of fp must be the number of beads
% Assume all the rois will have the same size:
roisize=fread(fp,2,clA,(nroi-1)*2);%Reads in the size of the first roi and then skips down the file to the data

grid_num = ceil(sqrt(nroi));

alldatatemp=fread(fp,inf,clA);%Reads in the rest of the file, which is 120 frames worth of intensities for all the rois
%First reshape into an x-by-y-by-nroi-by-fpw matrix:
alldata=reshape(alldatatemp,roisize(1),roisize(2),nroi,fpw);

for f=1:fpw
    clear currframetemp currframe
    currframetemp=reshape(alldata(:,:,:,f),size(alldata,1),size(alldata,2)*size(alldata,3)); %This puts all the beads for one frame into a long row;
    currframetemp=[currframetemp zeros(size(alldata,1),(grid_num^2-nroi)*size(alldata,2))];%zeros(0) returns an empty matrix, which is what we want; this pads the matrix with zeros for any empty grid spots
    %Can I avoid a for-loop here?
    for i=1:grid_num
        currframe(((i-1)*roisize(1)+1):i*roisize(1),:)=currframetemp(:,((i-1)*grid_num*roisize(2)+1):i*grid_num*roisize(2));
    end
    d{f} = [currframe];
    imagesc(currframe)
    colormap gray %Necessary?
    hold on
    plot([350 360], [350 350], '-w', 'LineWidth', 5)
    hold off
    text(352, 347, '1 micron', 'HorizontalAlignment', 'center')
    drawnow %Necessary?
end

fclose(fp);