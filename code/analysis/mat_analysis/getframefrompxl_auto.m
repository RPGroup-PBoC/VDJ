%frame=getframefrompxl(path,filenum,framenum,bdnum)
%
%Steph 10/14/08
%
%Same as getframfrompxl but user inputs the path to the directory with the
%pxl files instead of using uigetdir.
%
%Function to extract the desired frame from a .pxl TPM data file.  Based
%largely on open_pxl.m; from that code, the .pxl files are binary files
%where the first element is the number of frames per file (e.g. 120 for the
%Andor, 240 for the Basler), the second element is the number of beads
%(nroi), the third and fourth elements are the width and height of the roi
%for bead 1, fifth and sixth are for bead 1, etc., then there are
%width*height elements for bead 1, frame 1, then for bead 1, frame 1, etc.,
%then bead 1 frame 2, all the way through bead r, frame fpw.

function frame=getframefrompxl_auto(folder,filenum,framenum,bdnum)

if ismac
	ind=findstr(folder,'/');
    path=strcat(folder,'/');
else
	ind=findstr(folder,'\');
    path=strcat(folder,'\');
end
filename=folder(ind(end)+1:end);

fp = fopen(strcat(path,filename,'_',int2str(filenum),'.pxl'),'r');%Open file read-only
fpw = fread(fp,1,'uint16'); %Reads one element from the file fp, sets pointer at the second; first element must be ... ? number of frames per file I think
if framenum>fpw
    'Frame exceeds file dimensions.'
    return
end
nroi = fread(fp,1,'uint16');%Second element of fp must be the number of beads
if bdnum>nroi
    'Bead does not exist.'
    return
end
for r = 1:nroi
    roisize(r,1:2) = fread(fp,2,'uint16');%Elements three and four are "roisize" for bead one, five and six for bead two, etc
end
%Iterate through fp until reach the first element of framenum for bdnum
if framenum>1
    for f=1:framenum-1
        for r=1:nroi
            trash=fread(fp,roisize(r,1)*roisize(r,2),'uint16');
        end
    end    
end
%Now the pointer is at the first element for bd 1, framenum
if bdnum>1
    for r=1:bdnum-1
        trash=fread(fp,roisize(r,1)*roisize(r,2),'uint16');
    end
end
%Now the pointer is at the first element of bdnum,framenum
frame = fread(fp,roisize(r,1)*roisize(r,2),'uint16');%Reads roisize(1)*roisize(2) elements for bead r
frame = reshape(frame,roisize(r,1),roisize(r,2));%Now frame is an roisize(1)*roisize(2) matrix
figure,imagesc(frame)
colormap gray

fclose(fp);