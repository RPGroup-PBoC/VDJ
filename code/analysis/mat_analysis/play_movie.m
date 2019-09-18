%Steph 12/09

%Companion code to PXLtomatrix--displays the images loaded in PXLtomatrix
%as a 4-second movie.  Input is a matrix that is the output of PXLtomatrix.

function play_movie(alldata,roisize1,roisize2,nroi,fpw)

roisize(1)=roisize1;
roisize(2)=roisize2;

grid_num = ceil(sqrt(nroi));

figure

for f=1:fpw
    clear currframetemp currframe
    currframetemp=reshape(alldata(:,:,:,f),size(alldata,1),size(alldata,2)*size(alldata,3)); %This puts all the beads for one frame into a long row;
    currframetemp=[currframetemp zeros(size(alldata,1),(grid_num^2-nroi)*size(alldata,2))];%zeros(0) returns an empty matrix, which is what we want; this pads the matrix with zeros for any empty grid spots
    %Can I avoid a for-loop here?
    for i=1:grid_num
        currframe(((i-1)*roisize(1)+1):i*roisize(1),:)=currframetemp(:,((i-1)*grid_num*roisize(2)+1):i*grid_num*roisize(2));
    end
    imagesc(currframe)
    colormap gray %Necessary?
    drawnow %Necessary?
end
close