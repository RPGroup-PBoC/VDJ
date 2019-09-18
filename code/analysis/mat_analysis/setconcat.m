%Steph 12/5/08

%Scope has been drifting a lot, so have had to take more than two sets on
%some areas; this script concatenates sets by copying them into a new
%folder with file names that match that new folder; it also creates a new
%ANAL.mat file

%Preliminaries:
%concatfoldername: full path and name of new folder where the concatenated
%   sets will be stored--YOU MUST MAKE THIS FOLDER FIRST--do you? yes
%newname: first part of the new filenames to make; should be the last part
%   of concatfoldername
%foldtoconcat1, foldtoconcat2: full paths and names of folders to
%   concatenated; name1 and name2 are the same as newname but for the existing
%   folders
%CHANGE SAVENAME AT THE END

%concatfoldername='/Volumes/drobo/stephj/101001_TA94_oldscope_500fMSJlac_area1_1and2concat';
concatfoldername='T:\stephj\110506_TA82wprom_oldscope_100pMSJlac_area2_1and2concat';
newname='110506_TA82wprom_oldscope_100pMSJlac_area2_1and2concat';

%foldtoconcat1 = '/Volumes/drobo/stephj/101001_TA94_oldscope_500fMSJlac_area1';
foldtoconcat1 = 'T:\stephj\110506_TA82wprom_oldscope_100pMSJlac_area2_discprev';
name1 = '110506_TA82wprom_oldscope_100pMSJlac_area2_discprev';

%foldtoconcat2 = '/Volumes/drobo/stephj/101001_TA94_oldscope_500fMSJlac_area1_2';
foldtoconcat2 = 'T:\stephj\110506_TA82wprom_oldscope_100pMSJlac_area2_2';
name2='110506_TA82wprom_oldscope_100pMSJlac_area2_2';

%The code:

D1pxl=dir([foldtoconcat1,filesep,'*.pxl']);
D1mat=dir([foldtoconcat1,filesep,'*_POS.mat']);
anal1=load(strcat(foldtoconcat1,'/',name1,'ANAL.mat'));
pos1=anal1.rms_pos;
%******For cutting out some of the first set
% cutafterfile=875;
% temp=anal1.rms_pos;
% pos1=temp(:,1:cutafterfile);
%*******

D2pxl=dir([foldtoconcat2,filesep,'*.pxl']);
D2mat=dir([foldtoconcat2,filesep,'*_POS.mat']);
anal2=load(strcat(foldtoconcat2,'/',name2,'ANAL.mat'));
pos2=anal2.rms_pos;

for i=1:length(D1pxl) %Should be same number of pxls and mats
%for i=1:cutafterfile %******For cutting out some of the first set
    copyfile(strcat(foldtoconcat1,'/',name1,'_',int2str(i),'.pxl'),...
        strcat(concatfoldername,'/',newname,'_',int2str(i),'.pxl'));
    copyfile(strcat(foldtoconcat1,'/',name1,'_',int2str(i),'_POS.mat'),...
        strcat(concatfoldername,'/',newname,'_',int2str(i),'_POS.mat'));
    strcat(int2str(i),'/',int2str(length(D1pxl)))
end

for i=1:length(D2pxl) 
    copyfile(strcat(foldtoconcat2,'/',name2,'_',int2str(i),'.pxl'),...
        strcat(concatfoldername,'/',newname,'_',int2str(i+length(D1pxl)),'.pxl'));
    copyfile(strcat(foldtoconcat2,'/',name2,'_',int2str(i),'_POS.mat'),...
        strcat(concatfoldername,'/',newname,'_',int2str(i+length(D1pxl)),'_POS.mat'));
%*******For cutting out some of the first set:
% copyfile(strcat(foldtoconcat2,'/',name2,'_',int2str(i),'.pxl'),...
%         strcat(concatfoldername,'/',newname,'_',int2str(i+cutafterfile),'.pxl'));
% copyfile(strcat(foldtoconcat2,'/',name2,'_',int2str(i),'_POS.mat'),...
%         strcat(concatfoldername,'/',newname,'_',int2str(i+cutafterfile),'_POS.mat'));
%********
     strcat(int2str(i),'/',int2str(length(D2pxl)))
end

%rms_pos is the only part of the ANAL file that is used to load data:
rms_pos=[pos1 pos2];
cd(concatfoldername)
save 110506_TA82wprom_oldscope_100pMSJlac_area2_1and2concatANAL rms_pos %save strcat(newname,'ANAL') doesn't work ...

