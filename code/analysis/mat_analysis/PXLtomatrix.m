function alldata=PXLtomatrix(path, filename, filenumber)
    clA='uint16';
    
    imgfile = fopen(fullfile(path,strcat(filename,'_',int2str(filenumber),'.pxl')),'r');%Open file read-only
    fpw = fread(imgfile,1,clA); %Reads one element from the file imgfile, sets pointer at the second; first element is number of frames per file
    nroi = fread(imgfile,1,clA);%Second element is the number of rois
    roisizes = fread(imgfile,2,clA,(nroi-1)*2);
        
    alldatatemp=fread(imgfile,inf,clA);
    alldata=reshape(alldatatemp,roisizes(1),roisizes(2),nroi,fpw);
    fclose(imgfile);
    
%     img=alldata(:,:,bdnum,imgnumber);
%     maxp=min(min(img));
%     minp=max(max(img));
%     
%     %threshold = ((1-frac)*minp + frac*maxp);
    
end
    
   