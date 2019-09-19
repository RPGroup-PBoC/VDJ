%POSfromPXL(path,filename, number,varargin)
%
%For creating a POS file from a pxl file post-analysis (based on the
%relevant parts of realtime and subfunctions).  Path is the full path to the 
%directory containing the POS and PXL files.
%Filename is root file name, i.e. no number attached.
%
%9/10 added an optional input that is a directory to save the POS files to 
%which is different from the directory of the pxl files; there can be a 
%second optional input for thresholding the raw images before cross
%correlation.
%Also automatically loads parameters from the info file instead
%of hardwiring them into the code. 
%
%NOTE if you're using this on data taken with the TETHER1/TETHER2 version
%of the code, you'll need to check the defaults in the if/else statements
%at the begining of the code.
%
%Steph 5/09


%[1 2 3 7 9 11 13 14 16 17 18 21]  

%[1 2 3 6 7 8 9 12 16 17 18]

function POSfromPXL(path,filename, number,varargin)

    pm_threshold=0.1; %Threshold for finding centroid of cross correlation

    %Load relevant parameters from the infofile, if it exists
    if exist(fullfile(path,strcat(filename,'INFO.mat')))
        infofile = load(fullfile(path,strcat(filename,'INFO.mat')));
        if isfield(infofile,'clA')
            clA = infofile.clA;
        else
            clA = 'uint16';
        end
        if isfield(infofile,'nmpix')
            nmpix = infofile.nmpix;
        else
            nmpix = 100;
        end
    else
        clA = 'uint16';
        nmpix = 100;
    end

    %Need the first file for the cross-correlation:
    templatefile=fopen(fullfile(path,strcat(filename,'_1','.pxl')),'r');
    fpw = fread(templatefile,1,clA); %Reads one element from the file imgfile, sets pointer at the second; first element is number of frames per file
    nroi = fread(templatefile,1,clA);%Second element is the number of rois
    roisizes=fread(templatefile,2,clA,(nroi-1)*2);
    rms_pos = zeros(1,nroi);
    pm_xy1 = zeros(fpw,nroi,2);
    templates=zeros(roisizes(1),roisizes(2),nroi);
    alldatatemplatestemp=fread(templatefile,inf,clA);
    alldatatemplates=reshape(alldatatemplatestemp,roisizes(1),roisizes(2),nroi,fpw);
    fclose(templatefile);
    clear alldatatemplatestemp
    
    %If user wants to threshold the images
    if ~isempty(varargin) && length(varargin)>1
       alldatatemplatestemp=alldatatemplates; 
       clear alldatatemplates
       thresh = varargin{2};
       alldatatemplates = alldatatemplatestemp.*(alldatatemplatestemp>thresh);
    end
    
    for a=1:nroi
        roiffttempl = fft2(alldatatemplates(:,:,a,1));
        roiffttempl(1,1) = 0;
        templates(:,:,a) = roiffttempl;
    end

    %The file to analyze:
    imgfile = fopen(fullfile(path,strcat(filename,'_',int2str(number),'.pxl')),'r');%Open file read-only
    fpw = fread(imgfile,1,clA); %Reads one element from the file imgfile, sets pointer at the second; first element is number of frames per file
    nroi = fread(imgfile,1,clA);%Second element is the number of rois
    
    roisizes=fread(imgfile,2,clA,(nroi-1)*2);
    rms_pos = zeros(1,nroi);
    pm_xy1 = zeros(fpw,nroi,2);
        
    alldatatemp=fread(imgfile,inf,clA);
    alldata=reshape(alldatatemp,roisizes(1),roisizes(2),nroi,fpw);
    fclose(imgfile);
    
    clear alldatatemp
    
    %If user wants to threshold the images
    if ~isempty(varargin) && length(varargin)>1
       alldatatemp=alldata; 
       clear alldata
       thresh = varargin{2};
       alldata = alldatatemp.*(alldatatemp>thresh);
    end
    
    for acq_fr=1:fpw 

        for i = 1:nroi
            clear roidata
            roidata=alldata(:,:,i,acq_fr);
            
            roifft = fft2(roidata);
            roifft(1,1) = 0;
            pm_xy1(acq_fr,i,:) = nmpix* calc_pm_xy1_MJ(roifft, templates(:,:,i), pm_threshold);
        end
    end

    mean_pos = reshape(mean(pm_xy1,1),nroi,2);
    for i = 1:nroi
        mot = squeeze(var(pm_xy1(:,i,:))); %computes variance for each dimension
        rms_pos(i) = sqrt(sum(mot));  % computes in-plane motion
    end
    %rms_pos
    pm_xy1(:,1,:);
    %pause
    
    if ~isempty(varargin)
        save(fullfile(varargin{1},strcat(filename,'_',int2str(number),'_POS')),'nroi','pm_xy1','mean_pos','rms_pos');
    else
        save(fullfile(path,strcat(filename,'_',int2str(number),'_POS')),'nroi','pm_xy1','mean_pos','rms_pos');
    end
end  