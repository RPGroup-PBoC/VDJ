%conctitr_SJlacI_CURRENT (script)
%
%Plot various concentration curves.  Updated from conctitr_SJlacI_better.m
%
%Steph 1/2011

%% Parameters for plotting
OpDist = 0; %if one, plot operator center-to-center distance on x-axs of 
    %length series; if 0, use loop length 

%% DEFAULT IN ALL CASES IS NO FIT
FitE8indiv = 0;
FitE8sub0indiv = 0;
FitE8sub02indiv = 0;

FitO1O1indiv = 0;
FitO1O1sub0indiv = 0;
FitO1O1sub02indiv = 0;

FitO2O1indiv = 0;
FitO2O1sub0indiv = 0;
FitO2O1sub02indiv = 0; %Use sub0 with sub02 for the others

FitE8glbl = 0; %This is E8OidO1O2 global fit
FitE8sub0glbl = 0;
FitE8sub02glbl = 0; %As of 3/2011 this uses sub0 for O2-O1

FitTAindiv = 0;
FitTAsub0indiv = 0;
FitTAsub02indiv = 0;

FitTAglbl = 0; %This is TA E8OidO1O2TA global fit
FitTAsub0glbl = 0; 
FitTAsub02glbl = 0; %As of 3/2011 this uses sub0 for O2-O1
FitTAsub02glbldimers = 0;

FitE107 = 0; %This is E107BvsM, enforcing E8OidO1O2TAsub02 Kd's

FitPUCtotal = 0;
FitPUCBvsM = 0; %global fit to the two looped states
FitPUCBvsMsameK = 0; %global fit with only one Kd
FitE8TAPUCBvsM = 0; %This is global fit to PUC B+M and all E8's and TA

FitE8TAPUC = 0; %global fit to all E8 data, TA, and PUC B vs M
FitE8TAPUCsub0 = 0;
FitE8TAPUCsub02 = 0;
FitE8TAPUCsub02dimers = 0; %This is the dimers model

FitOldindiv = 0; %Fit to the ``old'' lac batch
FitOldsub02indiv = 0;
FitNewindiv = 0; %Fit to the ``new'' lac batch
FitNewsub02indiv = 0;

% FITS: load fit parameters
if FitE8indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E894 conc curve/E8indivfit.mat');
    %[R,FitE8] = plotpLoopTheory(KidE8*10^-12,K1E8*10^-12,JE8*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfit.mat');
    [R,FitE8] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear KidE8 K1E8 JE8 SEKidE8 SEK1E8 SEJE8 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    FitE8indiv = 1;
end
if FitE8sub0indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E894 conc curve/E8indivfitsub0.mat');
    %[R,FitE8] = plotpLoopTheory(KidE8sub0*10^-12,K1E8sub0*10^-12,JE8sub0*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfitsub0.mat');
    [R,FitE8] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear KidE8sub0 K1E8sub0 JE8sub0 SEKidE8sub0 SEK1E8sub0 SEJE8sub0 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    FitE8indiv = 1;
end
if FitE8sub02indiv
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfitsub02.mat');
    [R,FitE8] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear KidE8sub0 K1E8sub0 JE8sub0 SEKidE8sub0 SEK1E8sub0 SEJE8sub0 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    FitE8indiv = 1;
end
if FitO1O1indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/O1O1indivfit.mat');
    %[R,FitO1O1] = plotpLoopTheory(K1O1O1*10^-12,K1O1O1*10^-12,JO1O1*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfit.mat');
    [R,FitO1O1] = plotpLoopTheory(K1*10^-12,K1*10^-12,J*10^-12);
    clear K1O1O1 JO1O1 SEK1O1O1 SEJO1O1 K1 J SEK1 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitO1O1sub0indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/O1O1indivfitsub0.mat');
    %[R,FitO1O1] = plotpLoopTheory(K1O1O1sub0*10^-12,K1O1O1sub0*10^-12,JO1O1sub0*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfitsub0.mat');
    [R,FitO1O1] = plotpLoopTheory(K1*10^-12,K1*10^-12,J*10^-12);
    clear K1O1O1sub0 JO1O1sub0 SEK1O1O1sub0 SEJO1O1sub0 K1 J SEK1 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitO1O1sub02indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/O1O1indivfitsub0.mat');
    %[R,FitO1O1] = plotpLoopTheory(K1O1O1sub0*10^-12,K1O1O1sub0*10^-12,JO1O1sub0*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfitsub02.mat');
    [R,FitO1O1] = plotpLoopTheory(K1*10^-12,K1*10^-12,J*10^-12);
    clear K1O1O1sub0 JO1O1sub0 SEK1O1O1sub0 SEJO1O1sub0 K1 J SEK1 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitO2O1indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/O2O1indivfit.mat');
    %[R,FitO2O1] = plotpLoopTheory(K1O2O1*10^-12,K2O2O1*10^-12,JO2O1*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfit.mat');
    [R,FitO2O1] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    clear K1O2O1 K2O2O1 JO2O1 SEK1O2O1 SEK2O2O1 SEJO2O1 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitO2O1sub0indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/O2O1indivfitsub0.mat');
    %[R,FitO2O1] = plotpLoopTheory(K1O2O1sub0*10^-12,K2O2O1sub0*10^-12,JO2O1sub0*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfitsub0.mat');
    [R,FitO2O1] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    clear K1O2O1sub0 K2O2O1sub0 JO2O1sub0 SEK1O2O1sub0 SEK2O2O1sub0 SEJO2O1sub0 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitO2O1sub02indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/O2O1indivfitsub0.mat');
    %[R,FitO2O1] = plotpLoopTheory(K1O2O1sub0*10^-12,K2O2O1sub0*10^-12,JO2O1sub0*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfitsub02.mat');
    [R,FitO2O1] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    clear K1O2O1sub0 K2O2O1sub0 JO2O1sub0 SEK1O2O1sub0 SEK2O2O1sub0 SEJO2O1sub0 K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    close
end
if FitE8glbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2fit.mat');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2fit.mat');
    JE8 = fitparams(1);
    Kid = fitparams(2);
    K1 = fitparams(3);
    K2 = fitparams(4);
    [R,FitE8global] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1global] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1global] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    clear Kid K1 K2 JE8 SEKid SEK1 SEK2 SEJE8 fitparams paramSEs paramnames paramsbs
end
if FitE8sub0glbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2sub0fit.mat');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2fitsub0.mat');
    JE8 = fitparams(1);
    Kid = fitparams(2);
    K1 = fitparams(3);
    K2 = fitparams(4);
    [R,FitE8global] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1global] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1global] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    clear Kid K1 K2 JE8 SEKid SEK1 SEK2 SEJE8 fitparams paramSEs paramnames paramsbs
    FitE8glbl=1;
end
if FitE8sub02glbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2sub0fit.mat');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2fitsub02.mat');
    JE8 = fitparams(1);
    Kid = fitparams(2);
    K1 = fitparams(3);
    K2 = fitparams(4);
    [R,FitE8global] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1global] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1global] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    clear Kid K1 K2 JE8 SEKid SEK1 SEK2 SEJE8 fitparams paramSEs paramnames paramsbs
    FitE8glbl=1;
end
if FitTAindiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/TAindivfitparams');
    %[R,FitTA] = plotpLoopTheory(KidTA*10^-12,K1TA*10^-12,JTA*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfit.mat');
    [R,FitTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear KidTA K1TA JTA SEKidTA SEK1TA SEJTA K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
end
if FitTAsub0indiv
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/TAsub0indivfitparams');
    %[R,FitTA] = plotpLoopTheory(Kid*10^-12,K1*10^-12,J*10^-12);
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfitsub0.mat');
    [R,FitTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear Kid K1 J SEKid SEK1 SEJ K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    FitTAindiv = 1;
end
if FitTAsub02indiv
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfitsub02.mat');
    [R,FitTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear Kid K1 J SEKid SEK1 SEJ K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
    FitTAindiv = 1;
end
if FitTAglbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2TAfit');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfit.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    Kid = fitparams(3);
    K1 = fitparams(4);
    K2 = fitparams(5);
    [R,FitTAglobal] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wTA] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wTA] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close 
    clear Kid K1 K2 JE8 JTA SEKid SEK1 SEK2 SEJE8 SEJTA fitparams paramSEs paramnames paramsbs
end
if FitTAsub0glbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2TAsub0fitparams');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfitsub0.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    Kid = fitparams(3);
    K1 = fitparams(4);
    K2 = fitparams(5);
    [R,FitTAglobal] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wTA] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wTA] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    clear Kid K1 K2 JE8 JTA SEKid SEK1 SEK2 SEJE8 SEJTA fitparams paramSEs paramnames paramsbs
    FitTAglbl=1;
end
if FitTAsub02glbl
    %load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2TAsub0fitparams');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    Kid = fitparams(3);
    K1 = fitparams(4);
    K2 = fitparams(5);
    [R,FitTAglobal] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wTA] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wTA] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wTA] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    clear Kid K1 K2 JE8 JTA SEKid SEK1 SEK2 SEJE8 SEJTA fitparams paramSEs paramnames paramsbs
    FitTAglbl=1;
end
if FitE107
    E107 = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/111014E8107BvsMfit.mat');
    Kds = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
    Kid = Kds.fitparams(3);
    K1 = Kds.fitparams(4);
    JB = E107.fitparams(1);
    JM = E107.fitparams(2);
    [R,FitE107B] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JB*10^-12,'twoloops',JM*10^-12);
    close
    [R,FitE107M] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JM*10^-12,'twoloops',JB*10^-12);
    close
    [R,FitE107tot] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JB+JM)*10^-12);
    close
    clear JB JM K1 Kid E107 Kds
end
if FitPUCBvsM
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/PUC306 conc curve/PUCBvsMfitparams')
    [R,FitPUCB] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JB*10^-12,'twoloops',JM*10^-12);
    close
    [R,FitPUCM] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JM*10^-12,'twoloops',JB*10^-12);
    close
    [R,FitPUCtot] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JB+JM)*10^-12);
    close
    clear JB JM K1 Kid SEJB SEJM SEK1 SEKid
end
if FitPUCBvsMsameK
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/PUC306 conc curve/PUCBvsMsameKfitparams')
    [R,FitPUCB] = plotpLoopTheory(Kid*10^-12,Kid*10^-12,JB*10^-12,'twoloops',JM*10^-12);
    close
    [R,FitPUCM] = plotpLoopTheory(Kid*10^-12,Kid*10^-12,JM*10^-12,'twoloops',JB*10^-12);
    close
    [R,FitPUCtot] = plotpLoopTheory(Kid*10^-12,Kid*10^-12,(JB+JM)*10^-12);
    close
    clear JB JM Kid SEJB SEJM SEKid
end
if FitPUCtotal
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/PUC306 conc curve/PUCtotfitparams')
    [R,FitPUC] = plotpLoopTheory(KidPUC*10^-12,K1PUC*10^-12,JPUC*10^-12);
    close
    clear KidPUC K1PUC JPUC SEKidPUC SEK1PUC SEJPUC
end
if FitE8TAPUC
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110218E8TAPUCfit.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    JPUCB = fitparams(3);
    JPUCM = fitparams(4);
    Kid = fitparams(5);
    K1 = fitparams(6);
    K2 = fitparams(7);
    [R,FitTAwPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wPUC] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wPUC] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    [R,FitPUCBwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCB*10^-12,'twoloops',JPUCM*10^-12);
    close
    [R,FitPUCMwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCM*10^-12,'twoloops',JPUCB*10^-12);
    close
    [R,FitPUCtotwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JPUCM+JPUCB)*10^-12);
    close
    clear JB JM Kid SEJB SEJM SEKid
    clear Kid K1 K2 JE8 JTA fitparams paramSEs paramnames paramsbs
end
if FitE8TAPUCsub0
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110218E8TAPUCfitsub0.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    JPUCB = fitparams(3);
    JPUCM = fitparams(4);
    Kid = fitparams(5);
    K1 = fitparams(6);
    K2 = fitparams(7);
    [R,FitTAwPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wPUC] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wPUC] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    [R,FitPUCBwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCB*10^-12,'twoloops',JPUCM*10^-12);
    close
    [R,FitPUCMwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCM*10^-12,'twoloops',JPUCB*10^-12);
    close
    [R,FitPUCtotwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JPUCM+JPUCB)*10^-12);
    close
    clear JB JM Kid SEJB SEJM SEKid
    clear Kid K1 K2 JE8 JTA fitparams paramSEs paramnames paramsbs
end
if FitE8TAPUCsub02
    %load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110218E8TAPUCfitsub02.mat');
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110324E8TAPUCsub02fit.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    JPUCB = fitparams(3);
    JPUCM = fitparams(4);
    Kid = fitparams(5);
    K1 = fitparams(6);
    K2 = fitparams(7);
    [R,FitTAwPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12);
    close
    [R,FitE8wPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO1O1wPUC] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12);
    close
    [R,FitO2O1wPUC] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12);
    close
    [R,FitPUCBwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCB*10^-12,'twoloops',JPUCM*10^-12);
    close
    [R,FitPUCMwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCM*10^-12,'twoloops',JPUCB*10^-12);
    close
    [R,FitPUCtotwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JPUCM+JPUCB)*10^-12);
    close
    clear JB JM Kid SEJB SEJM SEKid
    clear Kid K1 K2 JE8 JTA fitparams paramSEs paramnames paramsbs
end
if FitE8TAPUCsub02dimers
    load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110324E8TAPUCsub02dimersfit.mat');
    JE8 = fitparams(1);
    JTA = fitparams(2);
    JPUCB = fitparams(3);
    JPUCM = fitparams(4);
    Kid = fitparams(5);
    K1 = fitparams(6);
    K2 = fitparams(7);
    KDT = fitparams(8);
    [R,FitTAwPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JTA*10^-12,'dimers',KDT*10^-12);
    close
    [R,FitE8wPUC] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JE8*10^-12,'dimers',KDT*10^-12);
    close
    [R,FitO1O1wPUC] = plotpLoopTheory(K1*10^-12,K1*10^-12,JE8*10^-12,'dimers',KDT*10^-12);
    close
    [R,FitO2O1wPUC] = plotpLoopTheory(K1*10^-12,K2*10^-12,JE8*10^-12,'dimers',KDT*10^-12);
    close
    [R,FitPUCBwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCB*10^-12,'twoloops+dimers',JPUCM*10^-12,KDT*10^-12);
    close
    [R,FitPUCMwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,JPUCM*10^-12,'twoloops+dimers',JPUCB*10^-12,KDT*10^-12);
    close
    [R,FitPUCtotwET] = plotpLoopTheory(Kid*10^-12,K1*10^-12,(JPUCM+JPUCB)*10^-12,'dimers',KDT*10^-12);
    close
    clear JB JM Kid SEJB SEJM SEKid
    clear Kid K1 K2 JE8 JTA fitparams paramSEs paramnames paramsbs
    FitE8TAPUCsub02 = 1;
end
if FitOldindiv
    load('/Volumes/dumbo/stephj/TPM data analysis/Shortloops titrations BEST/110302OldLacindivfit.mat');
    [R,FitOld] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
end
if FitOldsub02indiv
    load('/Volumes/dumbo/stephj/TPM data analysis/Shortloops titrations BEST/110302OldLacsub02indivfit.mat');
    [R,FitOld] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
end
if FitNewindiv
    load('/Volumes/dumbo/stephj/TPM data analysis/Shortloops titrations BEST/110302NewLacindivfit.mat');
    [R,FitNew] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
end
if FitNewsub02indiv
    load('/Volumes/dumbo/stephj/TPM data analysis/Shortloops titrations BEST/110302NewLacsub02indivfit.mat');
    [R,FitNew] = plotpLoopTheory(K1*10^-12,K2*10^-12,J*10^-12);
    close
    clear K1 K2 J SEK1 SEK2 SEJ newstats newdistribs paramsbs totSEs totmeans
end

%% What data to plot:
JE8JTA = 0;
JE8JTAMeans3000 = 0;
JE8JTAMeanssub03000 = 0;
JE8JTAMeanssub023000 = 0;

E8allmeans = 0;
E8Means3000 = 0;
E8Meanssub03000 = 0;
E8Meanssub023000 = 0;

O1O1allmeans = 0;
O1O1Means3000 = 0;
O1O1Meanssub03000 = 0;
O1O1Meanssub023000 = 0;

O2O1allmeans = 0;
O2O1Means3000 = 0;
O2O1Meanssub03000 = 0; %As of 3/2011, decided to use O2-O1 sub0 with sub02 for the others
O2O1Meanssub023000 = 0; %Decided I don't like this

TAallmeans = 0;
TAMeans3000 = 0;
TAMeanssub03000 = 0;
TAMeanssub023000 = 0;

PUCallmeans = 0;
PUCMeans3000 = 0;

PUCBmeans = 0;
PUCMmeans = 0;
PUCBMeans3000 = 0;
PUCMMeans3000 = 0;

%Length series: assume if want E8, want TA as well
LengthsGauss = 0;
LengthsGaussB = 0; %Bottom state
LengthsGaussM = 0;
LengthsThresh = 0; %Threshholding method--as of 3/11 this is wrong!
LengthsThreshB = 0; 
LengthsThreshM = 0;
LengthsGauss3000 = 0; %Means3000, Gauss
LengthsGauss3000B = 0;
LengthsGauss3000M = 0;
LengthsThresh3000 = 0; %Means3000, Thresh 
LengthsThresh3000B = 0;
LengthsThresh3000M = 0;
LengthsThresh3000sub02 = 0;
LengthsThresh3000sub02B = 0;
LengthsThresh3000sub02M = 0; 

%Length series J-factors: these are relative to E894, and are based on
%thresholding, and more than 3000 seconds, only
LengthsJE8 = 0; 
LengthsJE8B = 0;
LengthsJE8M = 0;
LengthsJTA = 0; 
LengthsJTAB = 0;
LengthsJTAM = 0;
LengthsJE8sub02 = 0; 
LengthsJE8Bsub02 = 0;
LengthsJE8Msub02 = 0;
LengthsJTAsub02 = 0; 
LengthsJTABsub02 = 0;
LengthsJTAMsub02 = 0;
LengthsJE8BMratiosub02 = 0; %This plots JM/Jtot
LengthsJTABMratiosub02 = 0;

%E898 concentration curve, small bds
E898concsmallbdssub02Thresh = 0; %Right now only have the option for loading the thresholded, nonzeros subtracted, etc
E898concsmallbdssub02ThreshBM = 0;
bigbdstoo = 0; %This will also plot the 100 pM point with the bigger beads

%E8107 concentration curve, normal beads
E8107concsub02Thresh = 0;
E8107concsub02ThreshBM = 1;

%E8108 concentration curve, normal beads
E8108concsub02Thresh = 0;
E8108concsub02ThreshBM = 0;

%HGseqs, w prom:
HGwpromE83000 = 0;
HGwpromTA3000 = 0;
HGwpromE8sub023000 = 0;
HGwpromTAsub023000 = 0;
HGwpromE8sub023000M = 0;
HGwpromE8sub023000B = 0;
HGwpromTAsub023000M = 0;
HGwpromTAsub023000B = 0;

HGwpromE8sub023000J = 0;
HGwpromTAsub023000J = 0;
HGwpromE8sub023000JBM = 0;
HGwpromTAsub023000JBM = 0;
HGwpromE8sub023000JBMratio = 0;%This plots JM/Jtot
HGwpromTAsub023000JBMratio = 0;

HGwpromE8sub023000F = 0;
HGwpromTAsub023000F = 0;
HGwpromDeltaF = 0; %This plots FE8-FTA
HGwpromE8sub023000FBM = 0;
HGwpromTAsub023000FBM = 0;

%With and without promoter J's plotted together:
AlllengthsJssub023000 = 0;
AlllengthsJssub023000BvsM = 0; %Not implemented yet

%Controls
linkedchannels = 0;
LCMeans3000 = 0;
LCMeanssub03000 = 0;
LCMeanssub023000 = 0;

smallbds = 0;
smallbds3000 = 0;
smallbdssub03000 = 0;
smallbdssub023000 = 0;

%Other repressor batches
Lac2 = 0;
Lac23000 = 0;
Lac2sub03000 = 0;
Lac2sub023000 = 0;
OldLac = 0;
OldLac3000 = 0;
OldLacsub03000 = 0;
OldLacsub023000 = 0;
NewLac = 0;
NewLac3000 = 0;
NewLacsub03000 = 0;
NewLacsub023000 = 0;

%%
%%%%%%%%%%%%%% Below here is the code that runs based on the above inputs
%Load all the data needed

if JE8JTA || JE8JTAMeans3000 || JE8JTAMeanssub03000 || JE8JTAMeanssub023000
    %Old way to calculate Jratios:
%     if JE8JTA
%         [blank,dataE8] = LoadDataAnalysis('E894 conc curve','stats');
% 
%         pLoopsE8J = zeros(1,length(dataE8));
%         SEsE8J = zeros(1,length(dataE8));
% 
%         for i=1:length(dataE8)
%             pLoopsE8J(i)=dataE8{i}.pLoop;
%             SEsE8J(i)=dataE8{i}.SEpLoop;
%         end
%         
%         [blank,dataTA] = LoadDataAnalysis('TA94 conc curve','stats');
% 
%         pLoopsTAJ = zeros(1,length(dataTA));
%         SEsTAJ = zeros(1,length(dataTA));
% 
%         for i=1:length(dataTA)
%             pLoopsTAJ(i)=dataTA{i}.pLoop;
%             SEsTAJ(i)=dataTA{i}.SEpLoop;
%         end
%         
%         clear dataE8 dataTA
%         
%     elseif JE8JTAMeans3000
%         [blank,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','Means',3000);
%         pLoopsE8J = dataE8(:,1);
%         SEsE8J = dataE8(:,2);
%         
%         [blank,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','Means',3000);
%         pLoopsTAJ = dataTA(:,1);
%         SEsTAJ = dataTA(:,2);
%         
%         clear dataE8 dataTA
%         
%     elseif JE8JTAMeanssub03000
%         [blank,dataE8sub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
%         pLoopsE8J = dataE8sub0(:,1);
%         SEsE8J = dataE8sub0(:,2);
%         
%         [blank,dataTAsub0] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
%         pLoopsTAJ = dataTAsub0(:,1);
%         SEsTAJ = dataTAsub0(:,2);
%         
%         E8TAglblfitparams = load('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E8OidO1O2TAsub0fitparams.mat');
% 
%         avgJratio = E8TAglblfitparams.JTA/E8TAglblfitparams.JE8;
%         %Error on that:
%         avgJratioerr = sqrt((E8TAglblfitparams.SEJTA/E8TAglblfitparams.JE8)^2+...
%             (E8TAglblfitparams.JTA*E8TAglblfitparams.SEJE8/((E8TAglblfitparams.JE8)^2))^2);
% 
%         clear dataE8sub0 dataTAsub0 E8TAglblfitparams
%         
%     end
%     
%     concJTAE8 = 10^-12.*[0.5 1 2.5 5 10 50 100 500 10^3 10^4];
%     [JTAE8,SEsJTAE8] = calcJratio(pLoopsE8J([1 2 3 4 5 7 8 9 10 11]),pLoopsTAJ(3:end-1),...
%         SEsE8J([1 2 3 4 5 7 8 9 10 11]),SEsTAJ(3:end-1));
%         
%     JE8JTA = 1;
%     
% end

    %New way, with bootstrapped errors:
    if JE8JTA
        [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',1000);
        [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',1000);
    elseif JE8JTAMeans3000
        [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',3000);
        [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
        fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfit.mat');
    elseif JE8JTAMeanssub03000
        [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
        [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
        fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfitsub0.mat');
    elseif JE8JTAMeanssub023000
        [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
        [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
        fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
    end
    
    concJTAE8 = 10^-12.*[0.5 1 2.5 5 10 50 100 500 10^3 10^4];
    dataE8 = dataE8([1 2 3 4 5 7 8 9 10 11]);
    dataTA = dataTA(3:12);

    [dJ,errdJ,notneeded] = calcJratio_bootstrp(dataE8,dataTA);
    
    %Calculate the average Jratio and an error on that J from the global E8
    %and TA fit.  Note that JE8TA doesn't have a fit, but I don't use that
    %anyway.
    
    if ~JE8JTA
        avgJratio = fitparams.fitparams(2)/fitparams.fitparams(1);
        tempJratio = zeros(size(fitparams.paramsbs,1),1);
        for r = 1:size(fitparams.paramsbs,1)
            tempJratio(r) = fitparams.paramsbs(r,2)/fitparams.paramsbs(r,1);
        end
        avgJratioerr = std(tempJratio,1);
    end
    
    JE8JTA=1;
    
end

%E8 all means
if E8allmeans
    [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','stats');

    pLoopsE8 = zeros(1,length(dataE8));
    SEsE8 = zeros(1,length(dataE8));

    for i=1:length(dataE8)
        pLoopsE8(i)=dataE8{i}.pLoop;
        SEsE8(i)=dataE8{i}.SEpLoop;
    end
end

%E8 Means, >3000 seconds
if E8Means3000
    [concE8,dataE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','Means',3000);
    pLoopsE8 = dataE8(:,1);
    SEsE8 = dataE8(:,2);
end

%E8 Means Sub0, >3000 seconds
if E8Meanssub03000
    [concE8,dataE8sub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsE8sub0 = dataE8sub0(:,1);
    SEsE8sub0 = dataE8sub0(:,2);
end

%E8 Means Sub02, >3000 seconds
if E8Meanssub023000
    [concE8,dataE8sub02] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsE8sub0 = dataE8sub02(:,1);
    SEsE8sub0 = dataE8sub02(:,2);
end

%TA all means
if TAallmeans
    [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','stats');

    pLoopsTA = zeros(1,length(dataTA));
    SEsTA = zeros(1,length(dataTA));

    for i=1:length(dataTA)
        pLoopsTA(i)=dataTA{i}.pLoop;
        SEsTA(i)=dataTA{i}.SEpLoop;
    end
end

%TA Means, >3000 seconds
if TAMeans3000
    [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','Means',3000);
    pLoopsTA = dataTA(:,1);
    SEsTA = dataTA(:,2);
end


%TA Means Sub0, >3000 seconds
if TAMeanssub03000
    [concTA,dataTAsub0] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsTAsub0 = dataTAsub0(:,1);
    SEsTAsub0 = dataTAsub0(:,2);
end

%TA Means Sub02, >3000 seconds
if TAMeanssub023000
    [concTA,dataTAsub02] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsTAsub0 = dataTAsub02(:,1);
    SEsTAsub0 = dataTAsub02(:,2);
end

%O1E8O1 all means
if O1O1allmeans
    [concE8O1O1,dataE8O1O1] = LoadDataAnalysis('O1E894O1 conc curve','stats');

    pLoopsE8O1O1 = zeros(1,length(dataE8O1O1));
    SEsE8O1O1 = zeros(1,length(dataE8O1O1));

    for i=1:length(dataE8O1O1)
        pLoopsE8O1O1(i)=dataE8O1O1{i}.pLoop;
        SEsE8O1O1(i)=dataE8O1O1{i}.SEpLoop;
    end
end

%O1E8O1 Means, >3000 seconds
if O1O1Means3000
    [concE8O1O1,dataE8O1O1] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','Means',3000);
    pLoopsE8O1O1 = dataE8O1O1(:,1);
    SEsE8O1O1 = dataE8O1O1(:,2);
end

%O1E8O1 Means Sub0, >3000 seconds
if O1O1Meanssub03000
    [concE8O1O1,dataE8O1O1sub0] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsE8O1O1sub0 = dataE8O1O1sub0(:,1);
    SEsE8O1O1sub0 = dataE8O1O1sub0(:,2);
end

%O1E8O1 Means Sub02, >3000 seconds
if O1O1Meanssub023000
    [concE8O1O1,dataE8O1O1sub02] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsE8O1O1sub0 = dataE8O1O1sub02(:,1);
    SEsE8O1O1sub0 = dataE8O1O1sub02(:,2);
end

%O2E8O1 all means
if O2O1allmeans
    [concE8O2O1,dataE8O2O1] = LoadDataAnalysis('O2E894O1 conc curve','stats');

    pLoopsE8O2O1 = zeros(1,length(dataE8O2O1));
    SEsE8O2O1 = zeros(1,length(dataE8O2O1));

    for i=1:length(dataE8O2O1)
        pLoopsE8O2O1(i)=dataE8O2O1{i}.pLoop;
        SEsE8O2O1(i)=dataE8O2O1{i}.SEpLoop;
    end
end

%O2E8O1 Means, >3000 seconds
if O2O1Means3000
    [concE8O2O1,dataE8O2O1] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','Means',3000);
    pLoopsE8O2O1 = dataE8O2O1(:,1);
    SEsE8O2O1 = dataE8O2O1(:,2);
end

%O2E8O1 Means Sub0, >3000 seconds
if O2O1Meanssub03000
    [concE8O2O1,dataE8O2O1sub0] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsE8O2O1sub0 = dataE8O2O1sub0(:,1);
    SEsE8O2O1sub0 = dataE8O2O1sub0(:,2);
end

%O2E8O1 Means Sub02, >3000 seconds
if O2O1Meanssub023000
    [concE8O2O1,dataE8O2O1sub02] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsE8O2O1sub0 = dataE8O2O1sub02(:,1);
    SEsE8O2O1sub0 = dataE8O2O1sub02(:,2);
end

%PUC (total/bottom/middle) means
if PUCBmeans || PUCMmeans || PUCallmeans
    [concPUC,dataPUC] = LoadDataAnalysis('PUC306 conc curve','stats');

    pLoopsPUC = zeros(1,length(dataPUC));
    pLoopsPUCB = zeros(1,length(dataPUC));
    pLoopsPUCM = zeros(1,length(dataPUC));
    SEsPUC = zeros(1,length(dataPUC));
    SEsPUCB = zeros(1,length(dataPUC));
    SEsPUCM = zeros(1,length(dataPUC));

    for i=1:length(dataPUC)
        pLoopsPUC(i)=dataPUC{i}.pLoop;
        SEsPUC(i)=dataPUC{i}.SEpLoop;
        pLoopsPUCB(i)=dataPUC{i}.pB;
        SEsPUCB(i)=dataPUC{i}.SEpB;
        pLoopsPUCM(i)=dataPUC{i}.pM;
        SEsPUCM(i)=dataPUC{i}.SEpM;
    end
    clear dataPUC
end

%PUC total loop Means, >3000 seconds
if PUCMeans3000
    [concPUC,dataPUC] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','Means',3000);
    pLoopsPUC = dataPUC(:,1);
    SEsPUC = dataPUC(:,2);
end

%PUC bottom loop means >3000 seconds
if PUCBMeans3000
    [concPUC,dataPUCB] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsB',3000);
    pLoopsPUCB = zeros(1,length(dataPUCB));
    SEsPUCB = zeros(1,length(dataPUCB));
    for i=1:length(concPUC)
        pLoopsPUCB(i) = mean(dataPUCB{i});
        SEsPUCB(i) = std(dataPUCB{i})/sqrt(length(dataPUCB{i})-1);
    end
end
%PUC middle loop means >3000 seconds
if PUCMMeans3000
    [concPUC,dataPUCM] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsM',3000);
    pLoopsPUCM = zeros(1,length(dataPUCM));
    SEsPUCM = zeros(1,length(dataPUCM));
    for i=1:length(concPUC)
        pLoopsPUCM(i) = mean(dataPUCM{i});
        SEsPUCM(i) = std(dataPUCM{i})/sqrt(length(dataPUCM{i})-1);
    end
end

%Length series (total/bottom/middle) means, gaussian fitting method
if LengthsGauss || LengthsGaussB || LengthsGaussM
    [concLen,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','stats');

    pLoopsLenE8 = zeros(1,length(dataLen{1}));
    pLoopsLenTA = zeros(1,length(dataLen{2}));
    pLoopsLenE8B = zeros(1,length(dataLen{1}));
    pLoopsLenTAB = zeros(1,length(dataLen{2}));
    pLoopsLenE8M = zeros(1,length(dataLen{1}));
    pLoopsLenTAM = zeros(1,length(dataLen{2}));
    SEsLenE8 = zeros(1,length(dataLen{1}));
    SEsLenTA = zeros(1,length(dataLen{2}));
    SEsLenE8B = zeros(1,length(dataLen{1}));
    SEsLenTAB = zeros(1,length(dataLen{2}));
    SEsLenE8M = zeros(1,length(dataLen{1}));
    SEsLenTAM = zeros(1,length(dataLen{2}));
    
    for i=1:length(dataLen{1})
        pLoopsLenE8(i)=dataLen{1}{i}.pLoop;
        SEsLenE8(i)=dataLen{1}{i}.SEpLoop;
        pLoopsLenE8B(i)=dataLen{1}{i}.pB;
        SEsLenE8B(i)=dataLen{1}{i}.SEpB;
        pLoopsLenE8M(i)=dataLen{1}{i}.pM;
        SEsLenE8M(i)=dataLen{1}{i}.SEpM;
    end
    for i=1:length(dataLen{2})
        pLoopsLenTA(i)=dataLen{2}{i}.pLoop;
        SEsLenTA(i)=dataLen{2}{i}.SEpLoop;
        pLoopsLenTAB(i)=dataLen{2}{i}.pB;
        SEsLenTAB(i)=dataLen{2}{i}.SEpB;
        pLoopsLenTAM(i)=dataLen{2}{i}.pM;
        SEsLenTAM(i)=dataLen{2}{i}.SEpM;
    end
    clear dataLen
end

%Length series (total/bottom/middle) means, thresholding method
if LengthsThresh || LengthsThreshB || LengthsThreshM
    [concLen,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','Threshstats');

    pLoopsLenE8 = zeros(1,length(dataLen{1}));
    pLoopsLenTA = zeros(1,length(dataLen{2}));
    pLoopsLenE8B = zeros(1,length(dataLen{1}));
    pLoopsLenTAB = zeros(1,length(dataLen{2}));
    pLoopsLenE8M = zeros(1,length(dataLen{1}));
    pLoopsLenTAM = zeros(1,length(dataLen{2}));
    SEsLenE8 = zeros(1,length(dataLen{1}));
    SEsLenTA = zeros(1,length(dataLen{2}));
    SEsLenE8B = zeros(1,length(dataLen{1}));
    SEsLenTAB = zeros(1,length(dataLen{2}));
    SEsLenE8M = zeros(1,length(dataLen{1}));
    SEsLenTAM = zeros(1,length(dataLen{2}));
    
    for i=1:length(dataLen{1})
        pLoopsLenE8(i)=dataLen{1}{i}.pLoop;
        SEsLenE8(i)=dataLen{1}{i}.SEpLoop;
        pLoopsLenE8B(i)=dataLen{1}{i}.pB;
        SEsLenE8B(i)=dataLen{1}{i}.SEpB;
        pLoopsLenE8M(i)=dataLen{1}{i}.pM;
        SEsLenE8M(i)=dataLen{1}{i}.SEpM;
    end
    for i=1:length(dataLen{2})
        pLoopsLenTA(i)=dataLen{2}{i}.pLoop;
        SEsLenTA(i)=dataLen{2}{i}.SEpLoop;
        pLoopsLenTAB(i)=dataLen{2}{i}.pB;
        SEsLenTAB(i)=dataLen{2}{i}.SEpB;
        pLoopsLenTAM(i)=dataLen{2}{i}.pM;
        SEsLenTAM(i)=dataLen{2}{i}.SEpM;
    end
    clear dataLen
end

%Length series total loop Means, >3000 seconds, gaussian fitting method
if LengthsGauss3000
    [concLen,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','HistoAnal','SJLacI','Means',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8Gauss = dataLen{1}(:,1);
    pLoopsLenTAGauss = dataLen{2}(:,1);
    SEsLenE8Gauss = dataLen{1}(:,2);
    SEsLenTAGauss = dataLen{2}(:,2);
    clear concLen dataLen
end

%Length series total loop Means, >3000 seconds, thresholding method
if LengthsThresh3000
    [concLen,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Means',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8Thresh = dataLen{1}(:,1);
    pLoopsLenTAThresh = dataLen{2}(:,1);
    SEsLenE8Thresh = dataLen{1}(:,2);
    SEsLenTAThresh = dataLen{2}(:,2);
    clear concLen dataLen
end

%Length series total loop Means, >3000 seconds, thresholding method, sub02
if LengthsThresh3000sub02
    [concLen,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8Thresh = dataLen{1}(:,1);
    pLoopsLenTAThresh = dataLen{2}(:,1);
    SEsLenE8Thresh = dataLen{1}(:,2);
    SEsLenTAThresh = dataLen{2}(:,2);
    clear concLen dataLen
end

%Length series bottom and/or middle loop means >3000 seconds, gaussian fitting method
if LengthsGauss3000B || LengthsGauss3000M
    [concLen,dataLenB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','HistoAnal','SJLacI','MeansB',3000);
    [concLen,dataLenM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','HistoAnal','SJLacI','MeansM',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8GaussB = dataLenB{1}(:,1);
    pLoopsLenE8GaussM = dataLenM{1}(:,1);
    pLoopsLenTAGaussB = dataLenB{2}(:,1);
    pLoopsLenTAGaussM = dataLenM{2}(:,1);
    SEsLenE8GaussB = dataLenB{1}(:,2);
    SEsLenE8GaussM = dataLenM{1}(:,2);
    SEsLenTAGaussB = dataLenB{2}(:,2);
    SEsLenTAGaussM = dataLenM{2}(:,2);
    clear concLen dataLenB dataLenM
end

%Length series bottom and/or middle loop means >3000 seconds, thresholding method
if LengthsThresh3000B || LengthsThresh3000M
    [concLen,dataLenB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','MeansB',3000);
    [concLen,dataLenM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','MeansM',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8ThreshB = dataLenB{1}(:,1);
    pLoopsLenE8ThreshM = dataLenM{1}(:,1);
    pLoopsLenTAThreshB = dataLenB{2}(:,1);
    pLoopsLenTAThreshM = dataLenM{2}(:,1);
    SEsLenE8ThreshB = dataLenB{1}(:,2);
    SEsLenE8ThreshM = dataLenM{1}(:,2);
    SEsLenTAThreshB = dataLenB{2}(:,2);
    SEsLenTAThreshM = dataLenM{2}(:,2);
    clear concLen dataLenB dataLenM
end

%Length series bottom and/or middle loop means >3000 seconds, thresholding method, sub02
if LengthsThresh3000sub02B || LengthsThresh3000sub02M
    [concLen,dataLenB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02B',3000);
    [concLen,dataLenM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02M',3000);
    lengthsE8 = concLen{1};
    lengthsTA = concLen{2};
    pLoopsLenE8ThreshB = dataLenB{1}(:,1);
    pLoopsLenE8ThreshM = dataLenM{1}(:,1);
    pLoopsLenTAThreshB = dataLenB{2}(:,1);
    pLoopsLenTAThreshM = dataLenM{2}(:,1);
    SEsLenE8ThreshB = dataLenB{1}(:,2);
    SEsLenE8ThreshM = dataLenM{1}(:,2);
    SEsLenTAThreshB = dataLenB{2}(:,2);
    SEsLenTAThreshM = dataLenM{2}(:,2);
    clear concLen dataLenB dataLenM
end

%Length series J-factors
if LengthsJE8
    disp('FIX!')
    [Lens,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribs',3000);
    LensE8 = Lens{1};
    E8distrib = dataLen{1};
    clear Lens dataLen
    
    %Calculate relative to TA94, 100 pM, using the J-factor from the global
    %fit:
    fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfit.mat');
    JTA94 = fitparams.fitparams(2);
    clear fitparams
    
    %The error will be based on bootstrapping both the TA data and the
    %length series:
    [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    TAdistrib = dataTA{logical(concTA==100*10^-12)};
    clear concTA dataTA
    
    for i = 1:length(E8distrib)
        refdistrib{i} = TAdistrib; %Is there a faster way to do this?
    end
    clear TAdistrib
    
    [dJE8rel,errdJE8,notneeded] = calcJratio_bootstrp(E8distrib,refdistrib,JTA94);
    
end
if LengthsJTA
    disp('FIX!')
    [Lens,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribs',3000);
    LensTA = Lens{2};
    TAlendistrib = dataLen{2};
    clear Lens dataLen
    
    %Calculate relative to TA94, 100 pM, using the J-factor from the global
    %fit:
    fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110217E8OidO1O2TAfit.mat');
    JTA94 = fitparams.fitparams(2);
    clear fitparams
    
    %The error will be based on bootstrapping both the TA data and the
    %length series:
    [concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    TAdistrib = dataTA{logical(concTA==100*10^-12)};
    clear concTA dataTA
    
    for i = 1:length(TAlendistrib)
        refdistrib{i} = TAdistrib; %Is there a faster way to do this?
    end
    clear TAdistrib
    
    [dJTA,errdJTA,notneeded] = calcJratio_bootstrp(TAlendistrib,refdistrib,JTA94);

end

if LengthsJE8sub02 || LengthsJE8Bsub02 || LengthsJE8Msub02 || AlllengthsJssub023000 ||...
    LengthsJE8BMratiosub02 %Need the total J's for the separated pLoops
    [Lens,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribssub02',3000);
    LensE8 = Lens{1};
    E8distrib = dataLen{1};
    clear Lens dataLen
    
    fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
    
    %%%%% CHOOSE TO CALCULATE REL TO E894 OR TA94
    
    %*****Calculate relative to TA94, 100 pM, using the J-factor from the global
    %fit:
    %JTA94 = fitparams.fitparams(2);

    %*****Calculate relative to E894, 100 pM, using the J-factor from the global
    %fit:
    JTA94 = fitparams.fitparams(1)*10^-12; %Need this in M
%     ptot = getpLoopfromTheory(fitparams.fitparams(3).*10^-12,...
%         fitparams.fitparams(4).*10^-12,fitparams.fitparams(1).*10^-12,100
%         *10^-12); %Don't actually need this
    JTA94distrib = fitparams.paramsbs(:,1).*10^-12;
    
    %*****The error will be based on bootstrapping both the TA94 (or E894) data and the
    %length series:
    %[concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %[concTA,dataTA] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000); %calling it TA just so as not to have to rename variables
    %Update 3/2011: Calculate error relative to the 10^4 global fit results
    
    
    %%%%%
    
    TAdistrib = getpLoopfromTheory(fitparams.paramsbs(:,3).*10^-12,...
        fitparams.paramsbs(:,4).*10^-12,fitparams.paramsbs(:,1).*10^-12,100*10^-12);
    
    for i = 1:length(E8distrib)
        refdistribE{i} = TAdistrib; %Is there a faster way to do this?
    end
    clear TAdistrib
    
    [dJE8,errdJE8,notneeded] = calcJratio_bootstrp(E8distrib,refdistribE,JTA94,JTA94distrib);
    %3/2011: Decided to replace these relative J-factors for the 94 bp
    %lengths with the ones calculated from the fits
    dJE8(logical(LensE8==94)) = fitparams.fitparams(1)*10^-12;
    errdJE8(logical(LensE8==94)) = fitparams.paramSEs(1)*10^-12;
    clear fitparams E8distrib notneeded
    
end
if LengthsJTAsub02  || LengthsJTABsub02 || LengthsJTAMsub02|| AlllengthsJssub023000||...
    LengthsJTABMratiosub02
    [Lens,dataLen] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribssub02',3000);
    LensTA = Lens{2};
    TAlendistrib = dataLen{2};
    clear Lens dataLen
    
    fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
    
    %%%%% CHOOSE TO CALCULATE REL TO E894 OR TA94
    
    %*****Calculate relative to TA94, 100 pM, using the J-factor from the global
    %fit:
    %JTA94 = fitparams.fitparams(2);

    %*****Calculate relative to E894, 100 pM, using the J-factor from the global
    %fit:
    JTA94 = fitparams.fitparams(1)*10^-12;
%     ptot = getpLoopfromTheory(fitparams.fitparams(3).*10^-12,...
%         fitparams.fitparams(4).*10^-12,fitparams.fitparams(1).*10^-12,100
%         *10^-12); %Don't actually need this
    JTA94distrib = fitparams.paramsbs(:,1).*10^-12;
    
    %*****The error will be based on bootstrapping both the TA94 (or E894) data and the
    %length series:
    %[concTA,dataTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %update:see lengthsE8 section immediately previous 3.2011
    %[concTA,dataTA] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000); %calling it TA just so as not to have to rename variables
    
    %%%%%
    
    %TAdistrib = dataTA{logical(concTA==100*10^-12)};
    TAdistrib = getpLoopfromTheory(fitparams.paramsbs(:,3).*10^-12,...
        fitparams.paramsbs(:,4).*10^-12,fitparams.paramsbs(:,1).*10^-12,100*10^-12);
    
    for i = 1:length(TAlendistrib)
        refdistribT{i} = TAdistrib; %Is there a faster way to do this?
    end
    clear TAdistrib
    
    [dJTA,errdJTA,notneeded] = calcJratio_bootstrp(TAlendistrib,refdistribT,JTA94,JTA94distrib);
    %3/2011: Decided to replace these relative J-factors for the 94 bp
    %lengths with the ones calculated from the fits
    dJTA(logical(LensTA==94)) = fitparams.fitparams(2)*10^-12;
    errdJTA(logical(LensTA==94)) = fitparams.paramSEs(2)*10^-12;
    clear fitparams TAlendistrib notneeded
 
end

if LengthsJE8Bsub02 || LengthsJE8Msub02 || LengthsJTABsub02 || LengthsJTAMsub02 ||...
    LengthsJE8BMratiosub02 || LengthsJTABMratiosub02 %Assume if want one of B vs M, want the other; and also here assuming if want E8, want TA
    clear i
    oldrefdistrib = refdistribE;
    clear refdistribE
    refdistrib = oldrefdistrib{1};
    clear oldrefdistrib
    [notneeded,dataLenB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02B',3000);
    [notneeded,dataLenM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02M',3000);
    pLoopsLenE8B = dataLenB{1}(:,1);
    pLoopsLenE8M = dataLenM{1}(:,1);
    pLoopsLenTAB = dataLenB{2}(:,1);
    pLoopsLenTAM = dataLenM{2}(:,1);
    clear notneeded dataLenB dataLenM
    [notneeded,dataLenB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribssub02B',3000);
    [notneeded,dataLenM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','alldistribssub02M',3000);
    E8lendistribB = dataLenB{1};
    TAlendistribB = dataLenB{2};
    E8lendistribM = dataLenM{1};
    TAlendistribM = dataLenM{2};
    clear dataLenB dataLenM notneeded
    
    pBpME8 = pLoopsLenE8B./pLoopsLenE8M;
    pBpMTA = pLoopsLenTAB./pLoopsLenTAM;
    
    JME8 = dJE8./(pBpME8+1);
    JBE8 = dJE8 - JME8;
    JMTA = dJTA./(pBpMTA+1);
    JBTA = dJTA - JMTA;
    
    if LengthsJE8BMratiosub02
            %JBMratioE8 = JME8./JBE8;
            JBMratioE8 = JME8./dJE8;
    end
    if LengthsJTABMratiosub02
            %JBMratioTA = JMTA./JBTA;
            JBMratioTA = JMTA./dJTA;
    end
    
    %Bootstrap for the errors
    nboot = 10000;
    pBpME8bootstrp = zeros(nboot,length(LensE8));
    pBpMTAbootstrp = zeros(nboot,length(LensTA)); 
    JME8bstrp = zeros(nboot,length(LensE8));
    JBE8bstrp = zeros(nboot,length(LensE8));
    JMTAbstrp = zeros(nboot,length(LensTA));
    JBTAbstrp = zeros(nboot,length(LensTA));
    
    for i = 1:length(LensE8)
        tempvectorB = E8lendistribB{i};
        if size(tempvectorB,2) > 1
            tempvectorB = tempvectorB'; %data{i} needs to be a column vector
        end
        tempvectorM = E8lendistribM{i};
        if size(tempvectorM,2) > 1
            tempvectorM = tempvectorM'; %data{i} needs to be a column vector
        end
        temppLoopsB = zeros(nboot,length(tempvectorB));
        temppLoopsM = zeros(nboot,length(tempvectorB));
        %[temppLoopsB(:,i),newdistribs] = bootstrp(nboot,@(x)[mean(x)],tempvectorB); %Newdistribs will be a numbds x nboot matrix of resampling indices
        %[temppLoopsM(:,i)] = bootstrp(nboot,@(x)[mean(x)],tempvectorM); 
        %The above line is wrong!  Need to keep pB and pM
        %for each bead together!
        
        [notneeded,indices] = bootstrp(nboot,[],tempvectorB); %This just returns bootstrapped indices from 1 to length(tempvectorB)
        
        for ind = 1:nboot %is there a faster way?
            temppLoopsB(ind,:) = tempvectorB(indices(:,ind));
            temppLoopsM(ind,:) = tempvectorM(indices(:,ind));
        end
        
        temppLoopsBmeans = mean(temppLoopsB,2);
        temppLoopsMmeans = mean(temppLoopsM,2);
        
        %Switch back to the variable names used in the rest of this code:
        clear temppLoopsB temppLoopsM
        temppLoopsB(:,i) = temppLoopsBmeans;
        temppLoopsM(:,i) = temppLoopsMmeans;
        clear temppLoopsBmeans temppLoopsMmeans ind indices
        
        pBpME8bootstrp(:,i) = temppLoopsB(:,i)./temppLoopsM(:,i);
        ptotE8bootstrp(:,i) = temppLoopsB(:,i) + temppLoopsM(:,i);
        %Now I have a bootstrapped total looping probability; need a
        %bootstrapped Jtot, calculated relative to the refdistrib which is
        %from the global fits to the bootstrapped 94 bp data
        pratiotemp = (1-ptotE8bootstrp(:,i))./ptotE8bootstrp(:,i);
        pratioreftemp = (1-refdistrib)./refdistrib;
        dJE8bstrp(:,i) = JTA94distrib.*(pratioreftemp./pratiotemp);
        
        JME8bstrp(:,i) = dJE8bstrp(:,i)./(pBpME8bootstrp(:,i)+1);
        JBE8bstrp(:,i) = dJE8bstrp(:,i) - JME8bstrp(:,i);
        
        if LengthsJE8BMratiosub02
                %JBMratiobstrpE8(:,i) = JME8bstrp(:,i)./JBE8bstrp(:,i);
                JBMratiobstrpE8(:,i) = JME8bstrp(:,i)./dJE8bstrp(:,i);
        end
        
        clear tempvectorB tempvectorM temppLoopsB temppLoopsM pratiotemp pratioreftemp
    end
    clear i
    for i = 1:length(LensTA)
        tempvectorB = TAlendistribB{i};
        if size(tempvectorB,2) > 1
            tempvectorB = tempvectorB'; %data{i} needs to be a column vector
        end
        tempvectorM = TAlendistribM{i};
        if size(tempvectorM,2) > 1
            tempvectorM = tempvectorM'; %data{i} needs to be a column vector
        end
        temppLoopsB = zeros(nboot,length(tempvectorB));
        temppLoopsM = zeros(nboot,length(tempvectorB));
        %[temppLoopsB(:,i),newdistribs] = bootstrp(nboot,@(x)[mean(x)],tempvectorB); %Newdistribs will be a numbds x nboot matrix of resampling indices
        %[temppLoopsM(:,i)] = bootstrp(nboot,@(x)[mean(x)],tempvectorM); 
        %The above line is wrong!  Need to keep pB and pM
        %for each bead together!
        
        [notneeded,indices] = bootstrp(nboot,[],tempvectorB); %This just returns bootstrapped indices from 1 to length(tempvectorB)
        
        for ind = 1:nboot %is there a faster way?
            temppLoopsB(ind,:) = tempvectorB(indices(:,ind));
            temppLoopsM(ind,:) = tempvectorM(indices(:,ind));
        end
        
        temppLoopsBmeans = mean(temppLoopsB,2);
        temppLoopsMmeans = mean(temppLoopsM,2);
        
        %Switch back to the variable names used in the rest of this code:
        clear temppLoopsB temppLoopsM
        temppLoopsB(:,i) = temppLoopsBmeans;
        temppLoopsM(:,i) = temppLoopsMmeans;
        clear temppLoopsBmeans temppLoopsMmeans ind indices
        
        pBpMTAbootstrp(:,i) = temppLoopsB(:,i)./temppLoopsM(:,i);
        ptotTAbootstrp(:,i) = temppLoopsB(:,i) + temppLoopsM(:,i);
        %Now I have a bootstrapped total looping probability; need a
        %bootstrapped Jtot, calculated relative to the refdistrib which is
        %from the global fits to the bootstrapped 94 bp data
        pratiotemp = (1-ptotTAbootstrp(:,i))./ptotTAbootstrp(:,i);
        pratioreftemp = (1-refdistrib)./refdistrib;
        dJTAbstrp(:,i) = JTA94distrib.*(pratioreftemp./pratiotemp);
        
        JMTAbstrp(:,i) = dJTAbstrp(:,i)./(pBpMTAbootstrp(:,i)+1);
        JBTAbstrp(:,i) = dJTAbstrp(:,i) - JMTAbstrp(:,i);
        
        if LengthsJTABMratiosub02
                %JBMratiobstrpTA(:,i) = JMTAbstrp(:,i)./JBTAbstrp(:,i);
                JBMratiobstrpTA(:,i) = JMTAbstrp(:,i)./dJTAbstrp(:,i);
        end
        
        clear tempvectorB tempvectorM temppLoopsB temppLoopsM pratiotemp pratioreftemp
    end
    
    errJBE8 = std(JBE8bstrp,1,1);
    errJME8 = std(JME8bstrp,1,1);
    errJBTA = std(JBTAbstrp,1,1);
    errJMTA = std(JMTAbstrp,1,1);
    
    if LengthsJE8BMratiosub02
        errJBMratioE8 = std(JBMratiobstrpE8,1,1);
    end
    if LengthsJTABMratiosub02
        errJBMratioTA = std(JBMratiobstrpTA,1,1);
    end
    
    clear E8lendistribB E8lendistribM TAlendistribB TAlendistribM JBE8bstrp JBTAbstrp JME8bstrp JMTAbstrp JTA94
    clear JTA94distrib dJE8bstrp dJTAbstrp i nboot pBpME8 pBpME8bootstrp pBpMTA pBpMTAbootstrp pLoopsLenE8B
    clear pLoopsLenE8M pLoopsLenTAB pLoopsLenTAM ptotE8bootstrp ptotTAbootstrp refdistrib
    clear JBMratiobstrpE8 JBMratiobstrpTA
    
end

%E898 conc curve, small beads
if E898concsmallbdssub02Thresh
    [concE98,dataE98] = LoadDataAnalysis('E898conccurve','ThreshAnal','SJLacI','Meanssub0',3000);
    pLoopsE98 = dataE98(:,1);
    SEsE98 = dataE98(:,2);
    
    if bigbdstoo
        [lengths,data] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02',3000);
        allpLoopsbig = data{1}(:,1);
        pLoopBig = allpLoopsbig(logical(lengths{1}==98));
        allSEsbig = data{1}(:,2);
        SEBig = allSEsbig(logical(lengths{1}==98));
        clear lengths data allpLoopsbig allSEsbig
    end
    
end

%E898 conc curve, small beads, bottom and middle states separately
if E898concsmallbdssub02ThreshBM
    [concE98,dataE98B] = LoadDataAnalysis('E898conccurve','ThreshAnal','SJLacI','Meanssub0B',3000);
    [concE98,dataE98M] = LoadDataAnalysis('E898conccurve','ThreshAnal','SJLacI','Meanssub0M',3000);
    pLoopsE98B = dataE98B(:,1);
    SEsE98B = dataE98B(:,2);
    pLoopsE98M = dataE98M(:,1);
    SEsE98M = dataE98M(:,2);
    
    if bigbdstoo
        [lengths,dataB] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02B',3000);
        [notneeded,dataM] = LoadDataAnalysis('Oid-E,T 89 to 100-O1','ThreshAnal','SJLacI','Meanssub02M',3000);
        allpLoopsBbig = dataB{1}(:,1);
        pLoopBigB = allpLoopsBbig(logical(lengths{1}==98));
        allpLoopsMbig = dataM{1}(:,1);
        pLoopBigM = allpLoopsMbig(logical(lengths{1}==98));
        allSEsBbig = dataB{1}(:,2);
        SEBigB = allSEsBbig(logical(lengths{1}==98));
        allSEsMbig = dataM{1}(:,2);
        SEBigM = allSEsMbig(logical(lengths{1}==98));
        clear lengths dataB dataM allpLoopsBbig allpLoopsMbig allSEsBbig allSEsMbig
    end
    
end

%E8107 conc curve
if E8107concsub02Thresh
    [concE107,dataE107] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','Meanssub02',3000);
    pLoopsE107 = dataE107(:,1);
    SEsE107 = dataE107(:,2);

    clear dataE107
    
end

%E8107 conc curve, bottom and middle states separately
if E8107concsub02ThreshBM
    [concE107,dataE107B] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','Meanssub02B',3000);
    [concE107,dataE107M] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','Meanssub02M',3000);
    pLoopsE107B = dataE107B(:,1);
    SEsE107B = dataE107B(:,2);
    pLoopsE107M = dataE107M(:,1);
    SEsE107M = dataE107M(:,2);
    
    clear dataE107B dataE107M
    
end

%E8108 conc curve
if E8108concsub02Thresh
    [concE108,dataE108] = LoadDataAnalysis('E8108conccurve','ThreshAnal','SJLacI','Meanssub02',3000);
    pLoopsE108 = dataE108(:,1);
    SEsE108 = dataE108(:,2);

    clear dataE108
    
end

%E8108 conc curve, bottom and middle states separately
if E8108concsub02ThreshBM
    [concE108,dataE108B] = LoadDataAnalysis('E8108conccurve','ThreshAnal','SJLacI','Meanssub02B',3000);
    [concE108,dataE108M] = LoadDataAnalysis('E8108conccurve','ThreshAnal','SJLacI','Meanssub02M',3000);
    pLoopsE108B = dataE108B(:,1);
    SEsE108B = dataE108B(:,2);
    pLoopsE108M = dataE108M(:,1);
    SEsE108M = dataE108M(:,2);
    
    clear dataE108B
    
end

%HGseqs--being lazier than with other data sets about variable names ...
if HGwpromE83000 || HGwpromTA3000 || HGwpromE8sub023000 || HGwpromTAsub023000 || ...
        HGwpromE8sub023000J || HGwpromE8sub023000F || HGwpromDeltaF ||...
        AlllengthsJssub023000 || HGwpromE8sub023000JBM || HGwpromE8sub023000FBM ||...
        HGwpromE8sub023000JBMratio || HGwpromTAsub023000JBMratio %Need total for B vs M
    if HGwpromE83000 || HGwpromTA3000
        [conc,data] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','Means',3000);
    elseif HGwpromE8sub023000 || HGwpromTAsub023000 || HGwpromE8sub023000J ||...
            HGwpromE8sub023000F || HGwpromDeltaF || AlllengthsJssub023000 ||...
            HGwpromE8sub023000JBM || HGwpromE8sub023000FBM ||...
            HGwpromE8sub023000JBMratio || HGwpromTAsub023000JBMratio
        [conc,data] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','Meanssub02',3000);
    end
    lengthsE8 = conc{1};
    lengthsTA = conc{2};
    pLoopsE8 = data{1}(:,1);
    pLoopsTA = data{2}(:,1);
    SEsE8 = data{1}(:,2);
    SEsTA = data{2}(:,2);
    clear conc data
    
    if HGwpromE8sub023000J || HGwpromE8sub023000F || HGwpromDeltaF ||...
            AlllengthsJssub023000 || HGwpromE8sub023000JBM || HGwpromE8sub023000FBM||...
            HGwpromE8sub023000JBMratio || HGwpromTAsub023000JBMratio
        fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
        
        %Solve the pLoop equation for J, then use fitted Oid and O2 and 100
        %pM to find J's:
        Rwprom = 100*10^-12;
        ZminusJterm = 1+Rwprom/(fitparams.fitparams(3)*10^-12)+ ...
            Rwprom/(fitparams.fitparams(5)*10^-12)+...
            Rwprom^2/(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12);
        JE8wprom = pLoopsE8.*ZminusJterm./((1-pLoopsE8).*(0.5.*Rwprom./(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12)));
        JTAwprom = pLoopsTA.*ZminusJterm./((1-pLoopsTA).*(0.5.*Rwprom./(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12)));
    
        clear ZminusJterm
        
        %Errors: 
        %Fitparams contains bootstrapped fit parameters in paramsbs.
        %Bootstrap the looping probabilities:
        [conc,data] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','alldistribssub02',3000);
        distribE8 = data{1};
        distribTA = data{2};
        clear conc data
        
        nboot = 10^4;
        newmeansE8 = zeros(nboot,length(distribE8));
        newmeansTA = zeros(nboot,length(distribTA));
        %newdistribsE8 = cell(length(distribE8),1);
        %newdistribsTA = cell(length(distribTA),1);
        for k = 1:length(distribE8)
            tempvector = distribE8{k};
            if size(tempvector,2) > 1
                tempvector = tempvector'; %needs to be a column vector
            end
            %[newmeansE8(:,k),newdistribsE8{k}] = bootstrp(nboot,@(x)[mean(x)],tempvector); %Newdistribs{k} will be a numbds x nboot matrix of resampled distributions
            [newmeansE8(:,k),notneeded] = bootstrp(nboot,@(x)[mean(x)],tempvector); %Newdistribs{k} will be a numbds x nboot matrix of resampled distributions
            clear tempvector notneeded
        end
        clear k
        for kk = 1:length(distribTA)
            tempvector = distribTA{kk};
            if size(tempvector,2) > 1
                tempvector = tempvector'; %needs to be a column vector
            end
            %[newmeansTA(:,kk),newdistribsTA{kk}] = bootstrp(nboot,@(x)[mean(x)],tempvector); %Newdistribs{kk} will be a numbds x nboot matrix of resampled distributions
            [newmeansTA(:,kk),notneeded] = bootstrp(nboot,@(x)[mean(x)],tempvector); %Newdistribs{kk} will be a numbds x nboot matrix of resampled distributions
            clear tempvector notneeded
        end
        
        clear kk
        
        newZminusJterm = 1+Rwprom./(fitparams.paramsbs(:,3).*10^-12)+ ...
            Rwprom./(fitparams.paramsbs(:,5).*10^-12)+...
            Rwprom^2./(fitparams.paramsbs(:,3).*10^-12.*fitparams.paramsbs(:,5).*10^-12);
        for k = 1:nboot %Is there a way to avoid a for-loop?
            newJsE8(k,:) = newmeansE8(k,:).*newZminusJterm(k)./((1-newmeansE8(k,:)).*(0.5.*Rwprom./(fitparams.paramsbs(k,3).*10^-12.*fitparams.paramsbs(k,5).*10^-12)));
            newJsTA(k,:) = newmeansTA(k,:).*newZminusJterm(k)./((1-newmeansTA(k,:)).*(0.5.*Rwprom./(fitparams.paramsbs(k,3).*10^-12.*fitparams.paramsbs(k,5).*10^-12)));
        end
        clear k nboot
        errJE8wprom = std(newJsE8,1,1); %the newJsE8 and newJsTA should be nboot by lengths
        errJTAwprom = std(newJsTA,1,1);
        clear newZminusJterm newmeansE8 newmeansTA fitparams
        
        %To calculate Floop
        if HGwpromE8sub023000F || HGwpromDeltaF || HGwpromE8sub023000FBM
           
            FE8wprom = -log(JE8wprom);
            FTAwprom = -log(JTAwprom);
            
            %Errors:
            newFsE8 = -log(newJsE8);
            %For some reason some of the new F's for E8 are Inf, even though the rest 
                %are around 23ish.  Setting any infinite values to the avg (there are only
                %26 of them).
            tempbadFs = newFsE8(:,7);
            tempbadFs(find(isinf(tempbadFs))) = FE8wprom(7);
            newFsE8(:,7) = tempbadFs;
            
            newFsTA = -log(newJsTA);
            
            errFE8wprom = std(newFsE8,1,1);
            errFTAwprom = std(newFsTA,1,1);
            
            %To calculate DeltaF
            if HGwpromDeltaF
               
                FE8minusFTA = FE8wprom-FTAwprom; %There should be a TA for every E8 and vice versa
                
                %Errors:
                errFE8minusFTA = std(newFsE8-newFsTA,1,1);
                
            end
            
        end
        
    end
    
end

if HGwpromE8sub023000M || HGwpromE8sub023000B || HGwpromTAsub023000M || HGwpromTAsub023000B ||...
        HGwpromE8sub023000JBM || HGwpromE8sub023000FBM || HGwpromE8sub023000JBMratio || HGwpromTAsub023000JBMratio
    [conc,dataB] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','Meanssub02B',3000);
    [conc,dataM] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','Meanssub02M',3000);
    lengthsE8 = conc{1};
    lengthsTA = conc{2};
    pLoopsE8B = dataB{1}(:,1);
    pLoopsE8M = dataM{1}(:,1);
    pLoopsTAB = dataB{2}(:,1);
    pLoopsTAM = dataM{2}(:,1);
    SEsE8B = dataB{1}(:,2);
    SEsE8M = dataM{1}(:,2);
    SEsTAB = dataB{2}(:,2);
    SEsTAM = dataM{2}(:,2);
    clear conc dataB dataM
    
    if HGwpromE8sub023000JBM || HGwpromE8sub023000FBM || HGwpromE8sub023000JBMratio || HGwpromTAsub023000JBMratio
        clear i
       %Done similarly to the no-promoter case above
        pBpME8 = pLoopsE8B./pLoopsE8M;
        pBpMTA = pLoopsTAB./pLoopsTAM;
    
        JME8wprom = JE8wprom./(pBpME8+1);
        JBE8wprom = JE8wprom - JME8wprom;
        JMTAwprom = JTAwprom./(pBpMTA+1);
        JBTAwprom = JTAwprom - JMTAwprom;
        
        if HGwpromE8sub023000JBMratio
                %JBMratioE8wprom = JME8wprom./JBE8wprom;
                JBMratioE8wprom = JME8wprom./JE8wprom;
        end
        if HGwpromTAsub023000JBMratio
                %JBMratioTAwprom = JMTAwprom./JBTAwprom;
                JBMratioTAwprom = JMTAwprom./JTAwprom;
        end
    
        %Bootstrap for the errors
        nboot = 10000;
        pBpME8bootstrp = zeros(nboot,length(lengthsE8));
        pBpMTAbootstrp = zeros(nboot,length(lengthsTA)); 
        JME8bstrp = zeros(nboot,length(lengthsE8));
        JBE8bstrp = zeros(nboot,length(lengthsE8));
        JMTAbstrp = zeros(nboot,length(lengthsTA));
        JBTAbstrp = zeros(nboot,length(lengthsTA));
        ptotE8bootstrp = zeros(nboot,length(lengthsE8));
        ptotTAbootstrp = zeros(nboot,length(lengthsTA));
        dJE8bstrp = zeros(nboot,length(lengthsE8));
        dJTAbstrp = zeros(nboot,length(lengthsTA));
        JBMratiobstrpE8 = zeros(nboot,length(lengthsE8));
        JBMratiobstrpTA = zeros(nboot,length(lengthsTA));
        
        [notneeded,datadistribB] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','alldistribssub02B',3000);
        [notneeded,datadistribM] = LoadDataAnalysis('HGseqsWProm','ThreshAnal','SJLacI','alldistribssub02M',3000);
        distribE8B = datadistribB{1};
        distribTAB = datadistribB{2};
        distribE8M = datadistribM{1};
        distribTAM = datadistribM{2};
        clear notneeded datadistribB datadistribM
        
        fitparams = load('/Volumes/dumbo/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');
        Rwprom = 100*10^-12;
        ZminusJterm = 1+Rwprom/(fitparams.fitparams(3)*10^-12)+ ...
                Rwprom/(fitparams.fitparams(5)*10^-12)+...
                Rwprom^2/(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12);
    
        for i = 1:length(lengthsE8)
            tempvectorB = distribE8B{i};
            if size(tempvectorB,2) > 1
                tempvectorB = tempvectorB'; %data{i} needs to be a column vector
            end
            tempvectorM = distribE8M{i};
            if size(tempvectorM,2) > 1
                tempvectorM = tempvectorM'; %data{i} needs to be a column vector
            end
            temppLoopsB = zeros(nboot,length(tempvectorB));
            temppLoopsM = zeros(nboot,length(tempvectorM));

            [notneeded,indices] = bootstrp(nboot,[],tempvectorB); %This just returns bootstrapped indices from 1 to length(tempvectorB)

            for ind = 1:nboot %is there a faster way?
                temppLoopsB(ind,:) = tempvectorB(indices(:,ind));
                temppLoopsM(ind,:) = tempvectorM(indices(:,ind));
            end

            temppLoopsBmeans = mean(temppLoopsB,2);
            temppLoopsMmeans = mean(temppLoopsM,2);

            %Switch back to the variable names used in the rest of this code:
            clear temppLoopsB temppLoopsM
            temppLoopsB(:,i) = temppLoopsBmeans;
            temppLoopsM(:,i) = temppLoopsMmeans;
            clear temppLoopsBmeans temppLoopsMmeans ind indices

            pBpME8bootstrp(:,i) = temppLoopsB(:,i)./temppLoopsM(:,i);
            %For some reason some of these are NaN's, even though the rest 
                %are around 1.  Setting any NaN values to the avg (there aren't many that I've seen).
            tempbadPrs = pBpME8bootstrp(:,i);
            tempbadPrs(find(isnan(tempbadPrs))) = pBpME8(i);
            pBpME8bootstrp(:,i) = tempbadPrs;
            clear tempbadPrs
            %Apparently some can be inf's too:
            tempbadPrs = pBpME8bootstrp(:,i);
            tempbadPrs(find(isinf(tempbadPrs))) = pBpME8(i);
            pBpME8bootstrp(:,i) = tempbadPrs;
            clear tempbadPrs
            
            ptotE8bootstrp(:,i) = temppLoopsB(:,i) + temppLoopsM(:,i);
            %Now I have a bootstrapped total looping probability; need a
            %bootstrapped Jtot.  
            
            dJE8bstrp(:,i) = ptotE8bootstrp(:,i).*ZminusJterm./((1-ptotE8bootstrp(:,i)).*(0.5.*Rwprom./(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12)));
            
            JME8bstrp(:,i) = dJE8bstrp(:,i)./(pBpME8bootstrp(:,i)+1);
            JBE8bstrp(:,i) = dJE8bstrp(:,i) - JME8bstrp(:,i);
            clear tempvectorB tempvectorM temppLoopsB temppLoopsM ptotE8bootstrp
            
            if HGwpromE8sub023000JBMratio
                %JBMratiobstrpE8(:,i) = JME8bstrp(:,1)./JBE8bstrp(:,1);
                JBMratiobstrpE8(:,i) = JME8bstrp(:,i)./dJE8bstrp(:,i);
                %Because dJE8bstrp can be zero, this can generate Inf's and
                %Nan's.  So:
                tempbadrs = JBMratiobstrpE8(:,i);
                if ~isempty(find(isnan(tempbadrs))) ||...
                        ~isempty(find(isinf(tempbadrs))) 
                    %Find finds Nan's and Inf's as nonzero!
                    goodelems1 = tempbadrs(find(~isnan(tempbadrs)));
                    goodelems2 = tempbadrs(find(~isinf(goodelems1)));
                    tempbadrs(find(isnan(tempbadrs))) = mean(goodelems2);
                    tempbadrs(find(isinf(tempbadrs))) = mean(goodelems2);
                    JBMratiobstrpE8(:,i) = tempbadrs;
                    clear goodelems1 goodelems2
                end
                clear tempbadrs
            end
            
        end
        
        clear i
        
        for i = 1:length(lengthsTA)
            tempvectorB = distribTAB{i};
            if size(tempvectorB,2) > 1
                tempvectorB = tempvectorB'; %data{i} needs to be a column vector
            end
            tempvectorM = distribTAM{i};
            if size(tempvectorM,2) > 1
                tempvectorM = tempvectorM'; %data{i} needs to be a column vector
            end
            temppLoopsB = zeros(nboot,length(tempvectorB));
            temppLoopsM = zeros(nboot,length(tempvectorM));

            [notneeded,indices] = bootstrp(nboot,[],tempvectorB); %This just returns bootstrapped indices from 1 to length(tempvectorB)

            for ind = 1:nboot %is there a faster way?
                temppLoopsB(ind,:) = tempvectorB(indices(:,ind));
                temppLoopsM(ind,:) = tempvectorM(indices(:,ind));
            end

            temppLoopsBmeans = mean(temppLoopsB,2);
            temppLoopsMmeans = mean(temppLoopsM,2);

            %Switch back to the variable names used in the rest of this code:
            clear temppLoopsB temppLoopsM
            temppLoopsB(:,i) = temppLoopsBmeans;
            temppLoopsM(:,i) = temppLoopsMmeans;
            clear temppLoopsBmeans temppLoopsMmeans ind indices

            pBpMTAbootstrp(:,i) = temppLoopsB(:,i)./temppLoopsM(:,i);
            %For some reason some of these are NaN's, even though the rest 
                %are around 1.  Setting any NaN values to the avg (there aren't many that I've seen).
            tempbadPrs = pBpMTAbootstrp(:,i);
            tempbadPrs(find(isnan(tempbadPrs))) = pBpMTA(i);
            pBpMTAbootstrp(:,i) = tempbadPrs;
            clear tempbadPrs
            %Apparently some can be inf's too:
            tempbadPrs = pBpME8bootstrp(:,i);
            tempbadPrs(find(isinf(tempbadPrs))) = pBpME8(i);
            pBpME8bootstrp(:,i) = tempbadPrs;
            clear tempbadPrs
            
            ptotTAbootstrp(:,i) = temppLoopsB(:,i) + temppLoopsM(:,i);
            %Now I have a bootstrapped total looping probability; need a
            %bootstrapped Jtot
            
            dJTAbstrp(:,i) = ptotTAbootstrp(:,i).*ZminusJterm./((1-ptotTAbootstrp(:,i)).*(0.5.*Rwprom./(fitparams.fitparams(3)*10^-12*fitparams.fitparams(5)*10^-12)));
            
            JMTAbstrp(:,i) = dJTAbstrp(:,i)./(pBpMTAbootstrp(:,i)+1);
            JBTAbstrp(:,i) = dJTAbstrp(:,i) - JMTAbstrp(:,i);
            clear tempvectorB tempvectorM temppLoopsB temppLoopsM ptotTAbootstrp
            
            if HGwpromTAsub023000JBMratio
                %JBMratiobstrpTA(:,i) = JMTAbstrp(:,1)./JBTAbstrp(:,1);
                JBMratiobstrpTA(:,i) = JMTAbstrp(:,i)./dJTAbstrp(:,i);
                %Because dJE8bstrp can be zero, this can generate Inf's and
                %Nan's.  So:
                tempbadrs = JBMratiobstrpTA(:,i);
                %Find finds Nan's and Inf's as nonzero!
                goodelems1 = tempbadrs(find(~isnan(tempbadrs)));
                goodelems2 = tempbadrs(find(~isinf(goodelems1)));
                tempbadrs(find(isnan(tempbadrs))) = mean(goodelems2);
                tempbadrs(find(isinf(tempbadrs))) = mean(goodelems2);
                JBMratiobstrpTA(:,i) = tempbadrs;
                clear tempbadrs goodelems1 goodelems2
            end
            
        end
        
        errJBE8wprom = std(JBE8bstrp,1,1);
        errJME8wprom = std(JME8bstrp,1,1);
        errJBTAwprom = std(JBTAbstrp,1,1);
        errJMTAwprom = std(JMTAbstrp,1,1);
        
        if HGwpromE8sub023000JBMratio
            errJBMratioE8wprom = std(JBMratiobstrpE8,1,1);
        end
        if HGwpromTAsub023000JBMratio
            errJBMratioTAwprom = std(JBMratiobstrpTA,1,1);
        end
        
        if HGwpromE8sub023000FBM
            
            disp('Write code for Floop B vs M!')
            
        end

        clear distribE8B distribE8M distribTAB distribTAM JBE8bstrp JBTAbstrp JME8bstrp JMTAbstrp
        clear dJE8bstrp dJTAbstrp i nboot pBpME8 pBpME8bootstrp pBpMTA pBpMTAbootstrp pLoopsE8B
        clear pLoopsE8M pLoopsTAB pLoopsTAM ptotE8bootstrp ptotTAbootstrp ZminusJterm
        clear JBMratiobstrpE8 JBMratiobstrpTA

    end
end

%Linked channels, all means
if linkedchannels
    [concLC,dataLC] = LoadDataAnalysis('linked channels','stats');

    pLoopsLC = zeros(1,length(dataLC));
    SEsLC = zeros(1,length(dataLC));

    for i=1:length(dataLC)
        pLoopsLC(i)=dataLC{i}.pLoop;
        SEsLC(i)=dataLC{i}.SEpLoop;
    end
end

%Linked channels, Means, >3000 seconds
if LCMeans3000
    [concLC,dataLC] = LoadDataAnalysis('linked channels','HistoAnal','SJLacI','Means',3000);
    pLoopsLC = dataLC(:,1);
    SEsLC = dataLC(:,2);
end

%Linked channels, Means Sub0, >3000 seconds
if LCMeanssub03000
    [concLC,dataLCsub0] = LoadDataAnalysis('linked channels','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsLCsub0 = dataLCsub0(:,1);
    SEsLCsub0 = dataLCsub0(:,2);
end

%Linked channels, Means Sub02, >3000 seconds
if LCMeanssub023000
    [concLC,dataLCsub02] = LoadDataAnalysis('linked channels','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsLCsub0 = dataLCsub02(:,1);
    SEsLCsub0 = dataLCsub02(:,2);
end

%Old Lac all means
if smallbds
    [concSmBds,dataSmBds] = LoadDataAnalysis('Indicia Beads','stats');

    pLoopsSmBds = zeros(1,length(dataSmBds));
    SEsSmBds = zeros(1,length(dataSmBds));

    for i=1:length(dataSmBds)
        pLoopsSmBds(i)=dataSmBds{i}.pLoop;
        SEsSmBds(i)=dataSmBds{i}.SEpLoop;
    end
end

%Old lac Means, >3000 seconds
if smallbds3000
    [concSmBds,dataSmBds] = LoadDataAnalysis('Indicia Beads','HistoAnal','SJLacI','Means',3000);
    pLoopsSmBds = dataSmBds(:,1);
    SEsSmBds = dataSmBds(:,2);
end

%Small bds Means Sub0, >3000 seconds
if smallbdssub03000
    [concSmBds,dataSmBdssub0] = LoadDataAnalysis('Indicia Beads','HistoAnal','SJLacI','Meanssub0',3000);
    pLoopsSmBdssub0 = dataSmBdssub0(:,1);
    SEsSmBdssub0 = dataSmBdssub0(:,2);
end

%Small beads Means Sub02, >3000 seconds
if smallbdssub023000
    [concSmBds,dataSmBdssub02] = LoadDataAnalysis('Indicia Beads','HistoAnal','SJLacI','Meanssub02',3000);
    pLoopsSmBdssub0 = dataSmBdssub02(:,1);
    SEsSmBdssub0 = dataSmBdssub02(:,2);
end

%SJLac2 all means
if Lac2
    [concLac2,dataLac2] = LoadDataAnalysis('E894 conc curve','stats','SJLacI2');

    pLoopsLac2 = zeros(1,length(dataLac2));
    SEsLac2 = zeros(1,length(dataLac2));

    for i=1:length(dataLac2)
        pLoopsLac2(i)=dataLac2{i}.pLoop;
        SEsLac2(i)=dataLac2{i}.SEpLoop;
    end
end

%SJLac2 Means, >3000 seconds
if Lac23000
    [concLac2,dataLac2] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI2','Means',3000);
    pLoopsLac2 = dataLac2(:,1);
    SEsLac2 = dataLac2(:,2);
end

%SJLac2 Means Sub0, >3000 seconds
if Lac2sub03000
    [concLac2,dataLac2sub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI2','Meanssub0',3000);
    pLoopsLac2sub0 = dataLac2sub0(:,1);
    SEsLac2sub0 = dataLac2sub0(:,2);
end

%SJLac2 Means Sub02, >3000 seconds
if Lac2sub023000
    [concLac2,dataLac2sub02] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI2','Meanssub02',3000);
    pLoopsLac2sub0 = dataLac2sub02(:,1);
    SEsLac2sub0 = dataLac2sub02(:,2);
end

%New lac all means
if NewLac
    [concNew,dataNew] = LoadDataAnalysis('E894 conc curve','stats','New');

    pLoopsNew = zeros(1,length(dataNew));
    SEsNew = zeros(1,length(dataNew));

    for i=1:length(dataNew)
        pLoopsNew(i)=dataNew{i}.pLoop;
        SEsNew(i)=dataNew{i}.SEpLoop;
    end
end

%New lac Means, >3000 seconds
if NewLac3000
    [concNew,dataNew] = LoadDataAnalysis('E894 conc curve','HistoAnal','New','Means',3000);
    pLoopsNew = dataNew(:,1);
    SEsNew = dataNew(:,2);
end

%Newlac Means Sub0, >3000 seconds
if NewLacsub03000
    [concNew,dataNewsub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','New','Meanssub0',3000);
    pLoopsNewsub0 = dataNewsub0(:,1);
    SEsNewsub0 = dataNewsub0(:,2);
end

%New lac Means Sub02, >3000 seconds
if NewLacsub023000
    [concNew,dataNewsub02] = LoadDataAnalysis('E894 conc curve','HistoAnal','New','Meanssub02',3000);
    pLoopsNewsub0 = dataNewsub02(:,1);
    SEsNewsub0 = dataNewsub02(:,2);
end

%Old Lac all means
if OldLac
    [concOld,dataOld] = LoadDataAnalysis('E894 conc curve','stats','Old');

    pLoopsOld = zeros(1,length(dataOld));
    SEsOld = zeros(1,length(dataOld));

    for i=1:length(dataOld)
        pLoopsOld(i)=dataOld{i}.pLoop;
        SEsOld(i)=dataOld{i}.SEpLoop;
    end
end

%Old lac Means, >3000 seconds
if OldLac3000
    [concOld,dataOld] = LoadDataAnalysis('E894 conc curve','HistoAnal','Old','Means',3000);
    pLoopsOld = dataOld(:,1);
    SEsOld = dataOld(:,2);
end

%Old lac Means Sub0, >3000 seconds
if OldLacsub03000
    [concOld,dataOldsub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','Old','Meanssub0',3000);
    pLoopsOldsub0 = dataOldsub0(:,1);
    SEsOldsub0 = dataOldsub0(:,2);
end

%Old lac Means Sub02, >3000 seconds
if OldLacsub023000
    [concOld,dataOldsub02] = LoadDataAnalysis('E894 conc curve','HistoAnal','Old','Meanssub02',3000);
    pLoopsOldsub0 = dataOldsub02(:,1);
    SEsOldsub0 = dataOldsub02(:,2);
end


%%

%%%%%%%%%% Plots: I write this section as I need different plots, so it doesn't actually plot
%everything possible given certain configurations of the inputs
if JE8JTA
    figure
    PlotHandleJ=errorbarxyHG(concJTAE8,dJ,[],errdJ,[],[],'ok','k');
    if JE8JTAMeans3000 || JE8JTAMeanssub03000 || JE8JTAMeanssub023000
        hold on
        PlotHandleJ(end+1)=plot([10^-15 10^-5],[avgJratio avgJratio],'-k');
        PlotHandleJ(end+1)=plot([10^-15 10^-5],[avgJratio+avgJratioerr avgJratio+avgJratioerr],'--k');
        PlotHandleJ(end+1)=plot([10^-15 10^-5],[avgJratio-avgJratioerr avgJratio-avgJratioerr],'--k');
    end
    xlim([10^-15 10^-5])
    %ylim([0 35])
    set(gca,'XScale','log')
    xlabel('LacI Concentration (M)','FontSize',16)
    ylabel('J_{loop, TA}/J_{loop, E8}','FontSize',16)
    set(gca,'FontSize',14)
    set(gca,'XTick',[10^-15 10^-14 10^-13 10^-12 10^-11 10^-10 10^-9 10^-8 10^-7 10^-6 10^-5])

    StandardFigure(PlotHandleJ,gca)
elseif AlllengthsJssub023000
    figure
    PlotHandleJ=plot(LensE8,dJE8,'.k',LensTA,dJTA,'.r',lengthsE8+36,JE8wprom,'ok',lengthsTA+36,JTAwprom,'or');
    hold on
    PlotHandleJ(end+1:end+length(LensE8)+1)=errorbarxyHG(LensE8,dJE8,[],errdJE8,[],[],'.k','k');
    PlotHandleJ(end+1:end+length(LensTA)+1)=errorbarxyHG(LensTA,dJTA,[],errdJTA,[],[],'.r','r');
    PlotHandleJ(end+1:end+length(lengthsE8)+1)=errorbarxyHG(lengthsE8+36,JE8wprom,[],errJE8wprom,[],[],'ok','k');
    PlotHandleJ(end+1:end+length(lengthsTA)+1)=errorbarxyHG(lengthsTA+36,JTAwprom,[],errJTAwprom,[],[],'or','r');
    
    xlabel('Loop Length (bp)','FontSize',16)
    xlim([88 125])
    ylabel('J_{loop} (M)','FontSize',16)
    set(gca,'FontSize',14)
    set(gca,'YScale','log')
    legend('Oid-E8-O1','Oid-TA-O1','Oid-E8-(prom)-O2','Oid-TA-(prom)-O2')

    StandardFigure(PlotHandleJ,gca)
elseif (LengthsJE8 && LengthsJTA) || (LengthsJE8sub02 && LengthsJTAsub02)
    figure
    PlotHandleJ=plot(LensE8,dJE8,'.k');
    hold on
    PlotHandleJ(end+1)=plot(LensTA,dJTA,'.r');
    PlotHandleJ(end+1:end+length(LensE8)+1)=errorbarxyHG(LensE8,dJE8,[],errdJE8,[],[],'.k','k');
    PlotHandleJ(end+1:end+length(LensTA)+1)=errorbarxyHG(LensTA,dJTA,[],errdJTA,[],[],'.r','r');
    %If I want to plot the Cloutier+Widom data too:
    
    xlabel('Loop Length (bp)','FontSize',16)
    %xlim([109 121])
    xlim([88 117])
    ylabel('J_{loop} (M)','FontSize',16)
    set(gca,'FontSize',14)
    set(gca,'YScale','log')
    legend('Oid-E8-O1','Oid-TA-O1')

    StandardFigure(PlotHandleJ,gca)
elseif (LengthsJE8B && LengthsJE8M) || (LengthsJE8Bsub02 && LengthsJE8Msub02)
    figure
    PlotHandleJ=plot(LensE8,JME8,'x--k');
    hold on
    PlotHandleJ(end+1)=plot(LensE8,JBE8,'o:k');
    PlotHandleJ(end+1)=plot(LensTA,JMTA,'x--r');
    PlotHandleJ(end+1)=plot(LensTA,JBTA,'o:r');
    PlotHandleJ(end+1:end+length(LensE8)+1)=errorbarxyHG(LensE8,JME8,[],errJME8,[],[],'x--k','k');
    PlotHandleJ(end+1:end+length(LensE8)+1)=errorbarxyHG(LensE8,JBE8,[],errJBE8,[],[],'o:k','k');
    PlotHandleJ(end+1:end+length(LensTA)+1)=errorbarxyHG(LensTA,JMTA,[],errJMTA,[],[],'x--r','r');
    PlotHandleJ(end+1:end+length(LensTA)+1)=errorbarxyHG(LensTA,JBTA,[],errJBTA,[],[],'o:r','r');
    
    xlabel('Loop Length (bp)','FontSize',16)
    %xlim([109 121])
    xlim([88 125])
    ylabel('J_{loop} (M)','FontSize',16)
    set(gca,'FontSize',14)
    set(gca,'YScale','log')
    legend('Oid-E8-O1, M','Oid-E8-O1, B','Oid-TA-O1, M', 'Oid-TA-O1, B')

    StandardFigure(PlotHandleJ,gca)
    
elseif HGwpromE8sub023000JBM
    if OpDist
        xaxisE = lengthsE8+36+20.5;
        xaxisT = lengthsTA+36+20.5;
    else
        xaxisE = lengthsE8+36;
        xaxisT = lengthsTA+36;
    end
    figure
    PlotHandleJ=plot(xaxisE,JME8wprom,'x--k');
    hold on
    PlotHandleJ(end+1)=plot(xaxisE,JBE8wprom,'o:k');
    PlotHandleJ(end+1)=plot(xaxisT,JMTAwprom,'x--r');
    PlotHandleJ(end+1)=plot(xaxisT,JBTAwprom,'o:r');
    PlotHandleJ(end+1:end+length(lengthsE8)+1)=errorbarxyHG(xaxisE,JME8wprom,[],errJME8wprom,[],[],'x--k','k');
    PlotHandleJ(end+1:end+length(lengthsE8)+1)=errorbarxyHG(xaxisE,JBE8wprom,[],errJBE8wprom,[],[],'o:k','k');
    PlotHandleJ(end+1:end+length(lengthsTA)+1)=errorbarxyHG(xaxisT,JMTAwprom,[],errJMTAwprom,[],[],'x--r','r');
    PlotHandleJ(end+1:end+length(lengthsTA)+1)=errorbarxyHG(xaxisT,JBTAwprom,[],errJBTAwprom,[],[],'o:r','r');
    
    xlabel('Loop Length (bp)','FontSize',16)
    if OpDist
        xlim([55+36+20.5 89+36+20.5])
    else
        xlim([88 89+36])
    end
    ylabel('J_{loop} (M)','FontSize',16)
    set(gca,'FontSize',14)
    set(gca,'YScale','log')
    legend('Oid-E8-(prom)-O2, M','Oid-E8-(prom)-O2, B','Oid-TA-(prom)-O2, M', 'Oid-TA-(prom)-O2, B')

    StandardFigure(PlotHandleJ,gca)
    
elseif HGwpromE8sub023000JBMratio %This actually plots JM/Jtot for with and no promoter, 
        %in two figures, one for E8 and the other for TA
    if OpDist
        xaxisEnoprom = LensE8 + 20.5;
        xaxisTnoprom = LensTA + 20.5;
        xaxisEwprom = lengthsE8+36+20.5;
        xaxisTwprom = lengthsTA+36+20.5;
    else
        xaxisEnoprom = LensE8;
        xaxisTnoprom = LensTA;
        xaxisEwprom = lengthsE8+36;
        xaxisTwprom = lengthsTA+36;
    end
    figure
    PlotHandleJ=plot(xaxisEnoprom,JBMratioE8,'.-k');
    hold on
    PlotHandleJ(end+1)=plot(xaxisEwprom,JBMratioE8wprom,'o--k');
    PlotHandleJ(end+1)=plot([88 89+36+20.5],[1 1],'--b');
    PlotHandleJ(end+1)=plot([88 89+36+20.5],[0 0],'--b');
    PlotHandleJ(end+1:end+length(LensE8)+1)=errorbarxyHG(xaxisEnoprom,JBMratioE8,[],errJBMratioE8,[],[],'.-k','k');
    PlotHandleJ(end+1:end+length(lengthsE8)+1)=errorbarxyHG(xaxisEwprom,JBMratioE8wprom,[],errJBMratioE8wprom,[],[],'o--k','k');
    
    xlabel('Loop Length (bp)','FontSize',16)
    if OpDist
        xlim([55+36+20.5 89+36+20.5])
    else
        xlim([88 89+36])
    end
    ylim([-0.2 1.4])
    ylabel('J_{loop, M}/J_{loop, tot}','FontSize',16)
    set(gca,'FontSize',14)
    %set(gca,'YScale','log')
    legend('Oid-E8-O1','Oid-E8-(prom)-O2')

    StandardFigure(PlotHandleJ,gca)

    figure
    PlotHandleJ2=plot(xaxisTnoprom,JBMratioTA,'.-r');
    hold on
    PlotHandleJ2(end+1)=plot(xaxisTwprom,JBMratioTAwprom,'o--r');
    PlotHandleJ2(end+1)=plot([88 89+36+20.5],[1 1],'--b');
    PlotHandleJ2(end+1)=plot([88 89+36+20.5],[0 0],'--b');
    PlotHandleJ2(end+1:end+length(LensTA)+1)=errorbarxyHG(xaxisTnoprom,JBMratioTA,[],errJBMratioTA,[],[],'.-r','r');
    PlotHandleJ2(end+1:end+length(lengthsTA)+1)=errorbarxyHG(xaxisTwprom,JBMratioTAwprom,[],errJBMratioTAwprom,[],[],'o--r','r');
    
    xlabel('Loop Length (bp)','FontSize',16)
    if OpDist
        xlim([55+36+20.5 89+36+20.5])
    else
        xlim([88 89+36])
    end
    ylim([-0.2 1.4])
    ylabel('J_{loop, M}/J_{loop, tot}','FontSize',16)
    set(gca,'FontSize',14)
    %set(gca,'YScale','log')
    legend('Oid-TA-O1','Oid-TA-(prom)-O2')

    StandardFigure(PlotHandleJ2,gca)
    
end


%CONTROLS
if linkedchannels || LCMeans3000 %all data, or just >3000
    %Always plotted with equivalent E8 set and global fit (the fit with TA
    %and all operators)
    [hLC,PlotHandleLC] = plotconccurve({concE8,concLC}, {pLoopsE8,pLoopsLC},...
            {SEsE8,SEsLC},{'Oid-E894-O1','Oid-E894-O1, linked channels'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
elseif LCMeanssub03000 || LCMeanssub023000 
    [hLC,PlotHandleLC] = plotconccurve({concE8,concLC}, {pLoopsE8sub0,pLoopsLCsub0},...
            {SEsE8sub0,SEsLCsub0},{'Oid-E894-O1','Oid-E894-O1, linked channels'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
elseif smallbds || smallbds3000 %all data, or just >3000
    %Always plotted with equivalent E8 set and global fit (the fit with TA
    %and all operators)
    [hsmbds,PlotHandlesmbds] = plotconccurve({concE8,concSmBds}, {pLoopsE8,pLoopsSmBds},...
            {SEsE8,SEsSmBds},{'Oid-E894-O1, 490 nm diameter beads','Oid-E894-O1, 270 nm diameter beads'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
elseif smallbdssub03000 || smallbdssub023000
    [hsmbds,PlotHandlesmbds] = plotconccurve({concE8,concSmBds}, {pLoopsE8sub0,pLoopsSmBdssub0},...
            {SEsE8sub0,SEsSmBdssub0},{'Oid-E894-O1, 490 nm diameter beads','Oid-E894-O1, 270 nm diameter beads'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
%OTHER BATCHES
elseif Lac2 || Lac23000
    %Always plotted with equivalent SJLac set and global fit (the fit with TA
    %and all operators)
    [hLac2,PlotHandleLac2] = plotconccurve({concE8,concLac2}, {pLoopsE8,pLoopsLac2},...
            {SEsE8,SEsLac2},{'Oid-E894-O1, SJLacI','Oid-E894-O1, SJLacI2'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
elseif Lac2sub03000 || Lac2sub023000
    [hLac2,PlotHandleLac2] = plotconccurve({concE8,concLac2}, {pLoopsE8sub0,pLoopsLac2sub0},...
            {SEsE8sub0,SEsLac2sub0},{'Oid-E894-O1, SJLacI','Oid-E894-O1, SJLacI2'},...
            {R},{FitE8wTA},{'E8','om'},{'E8'});
elseif OldLac || OldLac3000 %Right now if have Old will also have New
    %Plotting this with E8's *individual* fit
    [hOld,PlotHandleOld] = plotconccurve({concE8,concOld,concNew}, {pLoopsE8,pLoopsOld,pLoopsNew},...
            {SEsE8,SEsOld,SEsNew},{'Oid-E894-O1, SJLacI','Oid-E894-O1, KMLacI','Oid-E894-O1, KMLacI2'},...
            {R,R,R},{FitE8, FitOld, FitNew},{'E8','om','xm'},{'E8','--m',':m'});
elseif OldLacsub03000 || OldLacsub023000
    %Plotting this with E8's *individual* fit
    if FitE8sub02indiv
        [hOld,PlotHandleOld] = plotconccurve({concE8,concOld,concNew}, {pLoopsE8sub0,pLoopsOldsub0,pLoopsNewsub0},...
                {SEsE8sub0,SEsOldsub0,SEsNewsub0},{'Oid-E894-O1, SJLacI','Oid-E894-O1, KMLacI','Oid-E894-O1, KMLacI2'},...
                {R,R,R},{FitE8, FitOld, FitNew},{'E8','om','xm'},{'E8','--m',':m'});
    %Plotted with the global fit including TA
    elseif FitTAsub02glbl
        [hOld,PlotHandleOld] = plotconccurve({concE8,concOld,concNew}, {pLoopsE8sub0,pLoopsOldsub0,pLoopsNewsub0},...
                {SEsE8sub0,SEsOldsub0,SEsNewsub0},{'Oid-E894-O1, SJLacI','Oid-E894-O1, KMLacI','Oid-E894-O1, KMLacI2'},...
                {R,R,R},{FitE8wTA, FitOld, FitNew},{'E8','om','xb'},{'E8','--m',':b'});
    end
%REAL DATA
%Length series
elseif LengthsGauss3000 && LengthsThresh3000 %compare Gaussian fitting vs thresholding for the length series
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5,lengthsTA+20.5,lengthsE8+20.5+0.75,lengthsTA+20.5+0.75},... %Offset the two data sets
            {pLoopsLenE8Gauss,pLoopsLenTAGauss,pLoopsLenE8Thresh,pLoopsLenTAThresh},...
            {SEsLenE8Gauss,SEsLenTAGauss,SEsLenE8Thresh,SEsLenTAThresh},...
            {'E8, Gauss','TA, Gauss','E8, Thresh','TA, Thresh'},...
            [],[],{'E8','TA','ok','xr'},[],'Lengths','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8,lengthsTA,lengthsE8+0.75,lengthsTA+0.75},... %Offset the two data sets
            {pLoopsLenE8Gauss,pLoopsLenTAGauss,pLoopsLenE8Thresh,pLoopsLenTAThresh},...
            {SEsLenE8Gauss,SEsLenTAGauss,SEsLenE8Thresh,SEsLenTAThresh},...
            {'E8, Gauss','TA, Gauss','E8, Thresh','TA, Thresh'},...
            [],[],{'E8','TA','ok','xr'},[],'Lengths','LoopLen');
    end
elseif LengthsGauss3000B && LengthsThresh3000B %compare Gaussian fitting vs thresholding, separate looped states
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5,lengthsE+20.58,lengthsTA+20.5,lengthsTA+20.5,...
            lengthsE8+20.5+0.75,lengthsE8+20.5+0.75,lengthsTA+20.5+0.75,lengthsTA+20.5+0.75},... %Offset the two data sets
            {pLoopsLenE8GaussM,pLoopsLenE8GaussB,pLoopsLenTAGaussM,pLoopsLenTAGaussB,...
            pLoopsLenE8ThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshM,pLoopsLenTAThreshB},...
            {SEsLenE8GaussM,SEsLenE8GaussB,SEsLenTAGaussM,SEsLenTAGaussB,...
            SEsLenE8ThreshM,SEsLenE8ThreshB,SEsLenTAThreshM,SEsLenTAThreshB},...
            {'E8, Gauss, M','E8, Gauss, B','TA, Gauss, M','TA, Gauss, B','E8, Thresh, M','E8, Thresh, B','TA, Thresh, M','TA, Thresh, B'},...
            [],[],{'xk','ok','xr','or','xk','ok','xr','or'},[],'Lengths','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8,lengthsE8,lengthsTA,lengthsTA,...
            lengthsE8+0.75,lengthsE8+0.75,lengthsTA+0.75,lengthsTA+0.75},... %Offset the two data sets
            {pLoopsLenE8GaussM,pLoopsLenE8GaussB,pLoopsLenTAGaussM,pLoopsLenTAGaussB,...
            pLoopsLenE8ThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshM,pLoopsLenTAThreshB},...
            {SEsLenE8GaussM,SEsLenE8GaussB,SEsLenTAGaussM,SEsLenTAGaussB,...
            SEsLenE8ThreshM,SEsLenE8ThreshB,SEsLenTAThreshM,SEsLenTAThreshB},...
            {'E8, Gauss, M','E8, Gauss, B','TA, Gauss, M','TA, Gauss, B','E8, Thresh, M','E8, Thresh, B','TA, Thresh, M','TA, Thresh, B'},...
            [],[],{'xk','ok','xr','or','xk','ok','xr','or'},[],'Lengths','LoopLen');
    end
elseif (LengthsThresh3000 || LengthsThresh3000sub02) && ...
        ~(LengthsThresh3000B || LengthsThresh3000sub02B) %total looping probabilities alone
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5,lengthsTA+20.5},...
            {pLoopsLenE8Thresh,pLoopsLenTAThresh},...
            {SEsLenE8Thresh,SEsLenTAThresh},...
            {'Oid-E8-O1','Oid-TA-O1'},...
            [],[],{'.-k','.-r'},[],'Lengths','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8,lengthsTA},...
            {pLoopsLenE8Thresh,pLoopsLenTAThresh},...
            {SEsLenE8Thresh,SEsLenTAThresh},...
            {'Oid-E8-O1','Oid-TA-O1'},...
            [],[],{'.-k','.-r'},[],'Lengths','LoopLen');
    end
elseif ~(LengthsThresh3000 || LengthsThresh3000sub02) && ...
        (LengthsThresh3000B || LengthsThresh3000sub02B) %separate states alone; assume if want one looped state, want both
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5,lengthsE8+20.5,lengthsTA+20.5,lengthsTA+20.5},...
            {pLoopsLenE8ThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshM,pLoopsLenTAThreshB},...
            {SEsLenE8ThreshM,SEsLenE8ThreshB,SEsLenTAThreshM,SEsLenTAThreshB},...
            {'Oid-E8-O1, M','Oid-E8-O1, B','Oid-TA-O1, M','Oid-TA-O1, B'},...
            [],[],{'x--k','o:k','x--r','o:r'},[],'Lengths','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8,lengthsE8,lengthsTA,lengthsTA},...
            {pLoopsLenE8ThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshM,pLoopsLenTAThreshB},...
            {SEsLenE8ThreshM,SEsLenE8ThreshB,SEsLenTAThreshM,SEsLenTAThreshB},...
            {'Oid-E8-O1, M','Oid-E8-O1, B','Oid-TA-O1, M','Oid-TA-O1, B'},...
            [],[],{'x--k','o:k','x--r','o:r'},[],'Lengths','LoopLen');
    end
elseif (LengthsThresh3000 || LengthsThresh3000sub02) && ...
        (LengthsThresh3000B || LengthsThresh3000sub02B) %total pLoop and separate states; assume if want one looped state, want both
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5,lengthsE8+20.5,lengthsTA+20.5,lengthsTA+20.5,lengthsE8+20.5,lengthsTA+20.5},...
            {pLoopsLenE8Thresh,pLoopsLenTAThresh,pLoopsLenE8ThreshM,pLoopsLenTAThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshB},...
            {SEsLenE8Thresh,SEsLenTAThresh,SEsLenE8ThreshM,SEsLenTAThreshM,SEsLenE8ThreshB,SEsLenTAThreshB},...
            {'E8, All','TA, All','E8, M','TA, M','E8, B','TA, B'},...
            [],[],{'.-k','.-r','x--k','o--k','x--r','o--r'},[],'Lengths','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8,lengthsE8,lengthsTA,lengthsTA,lengthsE8,lengthsTA},...
            {pLoopsLenE8Thresh,pLoopsLenTAThresh,pLoopsLenE8ThreshM,pLoopsLenTAThreshM,pLoopsLenE8ThreshB,pLoopsLenTAThreshB},...
            {SEsLenE8Thresh,SEsLenTAThresh,SEsLenE8ThreshM,SEsLenTAThreshM,SEsLenE8ThreshB,SEsLenTAThreshB},...
            {'E8, All','TA, All','E8, M','TA, M','E8, B','TA, B'},...
            [],[],{'.-k','.-r','x--k','o--k','x--r','o--r'},[],'Lengths','OpDist');
    end
%HGseqs
elseif (HGwpromE83000 || HGwpromE8sub023000) && ...
        ~(HGwpromE8sub023000M || HGwpromE8sub023000B) %total looping probabilities alone; assume if want E8, want TA too
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5+36,lengthsTA+20.5+36},...
            {pLoopsE8,pLoopsTA},...
            {SEsE8,SEsTA},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'LengthsHG','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36,lengthsTA+36},...
            {pLoopsE8,pLoopsTA},...
            {SEsE8,SEsTA},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'LengthsHG','LoopLen');
    end
elseif ~(HGwpromE83000 || HGwpromE8sub023000) && ...
        (HGwpromE8sub023000M || HGwpromE8sub023000B) %separate states alone; assume if want one looped state, want both
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5+36,lengthsE8+20.5+36,lengthsTA+20.5+36,lengthsTA+20.5+36},...
            {pLoopsE8M,pLoopsE8B,pLoopsTAM,pLoopsTAB},...
            {SEsE8M,SEsE8B,SEsTAM,SEsTAB},...
            {'Oid-E8-(prom)-O2, M','Oid-E8-(prom)-O2, B','Oid-TA-(prom)-O2, M','Oid-TA-(prom)-O2, B'},...
            [],[],{'x--k','o:k','x--r','o:r'},[],'LengthsHG','OpDist');
    else
                [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36,lengthsE8+36,lengthsTA+36,lengthsTA+36},...
            {pLoopsE8M,pLoopsE8B,pLoopsTAM,pLoopsTAB},...
            {SEsE8M,SEsE8B,SEsTAM,SEsTAB},...
            {'Oid-E8-(prom)-O2, M','Oid-E8-(prom)-O2, B','Oid-TA-(prom)-O2, M','Oid-TA-(prom)-O2, B'},...
            [],[],{'x--k','o:k','x--r','o:r'},[],'LengthsHG','LoopLen');
    end
%J-factors, both states
elseif HGwpromE8sub023000J && ~HGwpromE8sub023000JBM %assume if want E8, want TA too
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+20.5+36,lengthsTA+20.5+36},...
            {JE8wprom,JTAwprom},...
            {errJE8wprom,errJTAwprom},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'JHG','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36,lengthsTA+36},...
            {JE8wprom,JTAwprom},...
            {errJE8wprom,errJTAwprom},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'JHG','LoopLen');

    end
%F_loop, both states
elseif HGwpromE8sub023000F && ~HGwpromE8sub023000FBM %assume if want E8, want TA too
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36+20.5,lengthsTA+36+20.5},...
            {FE8wprom,FTAwprom},...
            {errFE8wprom,errFTAwprom},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'FHG','OpDist');
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36,lengthsTA+36},...
            {FE8wprom,FTAwprom},...
            {errFE8wprom,errFTAwprom},...
            {'Oid-E8-(prom)-O2','Oid-TA-(prom)-O2'},...
            [],[],{'.-k','.-r'},[],'FHG','LoopLen');
    end
elseif HGwpromDeltaF
    if OpDist
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36+20.5,lengthsTA+36+20.5},...
            {FE8minusFTA,zeros(length(lengthsTA),1)},...
            {errFE8minusFTA,[]},...
            {'F_{E8}-F_{TA}'},...
            [],[],{'.-k','--k'},[],'FHG','OpDist');
        ylim([-3 3])
        ylabel('\Delta F_{loop,E8}-\Delta F_{loop,TA}')
    else
        [hLen,PlotHandleLen] = plotconccurve({lengthsE8+36,lengthsTA+36},...
            {FTAminusFE8,zeros(length(lengthsTA),1)},...
            {errFTAminusFE8,[]},...
            {'F_{E8}-F_{TA}'},...
            [],[],{'.-k','--k'},[],'FHG','LoopLen');
        ylim([-3 3])
        ylabel('\Delta F_{loop,E8}-\Delta F_{loop,TA}')
    end

%Concentration curves
%E8107 conc curve
elseif E8107concsub02Thresh && E8107concsub02ThreshBM
    if FitE107
        [hE107,PlotHandleE107] = plotconccurve({concE107,concE107,concE107},...
        {pLoopsE107,pLoopsE107M,pLoopsE107B}, {SEsE107,SEsE107M,SEsE107B},...
        {'Oid-E8107-O1,B+M','Oid-E8107-O1,M','Oid-E8107-O1,B'},...
        {R,R,R},{FitE107tot,FitE107M,FitE107B},{'.k','ok','xk'},{'-k','--k',':k'});
    else
        [hE107,PlotHandleE107] = plotconccurve({concE107,concE107,concE107},...
            {pLoopsE107,pLoopsE107M,pLoopsE107B}, {SEsE107,SEsE107M,SEsE107B},...
            {'Oid-E8107-O1,B+M','Oid-E8107-O1,M','Oid-E8107-O1,B'},...
            [],[],{'.k','ok','xk'});
    end
%E8108 conc curve
elseif E8108concsub02Thresh && E8108concsub02ThreshBM
    if FitE8108
        [hE108,PlotHandleE108] = plotconccurve({concE108,concE108,concE108},...
        {pLoopsE108,pLoopsE108M,pLoopsE108B}, {SEsE108,SEsE108M,SEsE108B},...
        {'Oid-E8108-O1,B+M','Oid-E8108-O1,M','Oid-E8108-O1,B'},...
        {R,R,R},{FitE108,FitE108M,FitE108B},{'.k','ok','xk'},{'-k','--k',':k'});
    else
        [hE108,PlotHandleE108] = plotconccurve({concE108,concE108,concE108},...
            {pLoopsE108,pLoopsE108M,pLoopsE108B}, {SEsE108,SEsE108M,SEsE108B},...
            {'Oid-E8108-O1,B+M','Oid-E8108-O1,M','Oid-E8108-O1,B'},...
            [],[],{'.k','ok','xk'});
    end

%E898 conc curve
elseif E898concsmallbdssub02Thresh && E898concsmallbdssub02ThreshBM
    if FitE898
        if bigbdstoo
            [hE98,PlotHandleE98] = plotconccurve({concE98,concE98,concE98,[100*10^-12],[100*10^-12],[100*10^-12]},...
                {pLoopsE98,pLoopsE98M,pLoopsE98B,pLoopBig,pLoopBigM,pLoopBigB}, {SEsE98,SEsE98M,SEsE98B,SEBig,SEBigM,SEBigB},...
                {'Oid-E898-O1,B+M','Oid-E898-O1,M','Oid-E898-O1,B'},...
                {R,R,R},{FitE98,FitE98M,FitE98B},{'.k','ok','xk','.m','om','xm'},{'-k',':k','--k'});
        else
            [hE98,PlotHandleE98] = plotconccurve({concE98,concE98,concE98},...
            {pLoopsE98,pLoopsE98M,pLoopsE98B}, {SEsE98,SEsE98M,SEsE98B},...
            {'Oid-E898-O1,B+M','Oid-E898-O1,M','Oid-E898-O1,B'},...
            {R,R,R},{FitE98,FitE98M,FitE98B},{'.k','ok','xk'},{'-k','--k',':k'});
        end
    else
        if bigbdstoo
            [hE98,PlotHandleE98] = plotconccurve({concE98,concE98,concE98,[100*10^-12],[100*10^-12],[100*10^-12]},...
                {pLoopsE98,pLoopsE98M,pLoopsE98B,pLoopBig,pLoopBigM,pLoopBigB}, {SEsE98,SEsE98M,SEsE98B,SEBig,SEBigM,SEBigB},...
                {'Oid-E898-O1,B+M','Oid-E898-O1,M','Oid-E898-O1,B'},...
                [],[],{'.k','ok','xk','.m','om','xm'});
        else
            [hE98,PlotHandleE98] = plotconccurve({concE98,concE98,concE98},...
                {pLoopsE98,pLoopsE98M,pLoopsE98B}, {SEsE98,SEsE98M,SEsE98B},...
                {'Oid-E898-O1,B+M','Oid-E898-O1,M','Oid-E898-O1,B'},...
                [],[],{'.k','ok','xk'});
        end
    end
    
elseif E8allmeans && ~O1O1allmeans && ~O2O1allmeans && ~TAallmeans && ~PUCMeans3000  %E8 alone, means
    if FitE8indiv || FitE8glbl
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, {R},{FitE8sub0},{'E8'},{'E8'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, [],[],{'E8'});
    end
elseif E8Means3000 && ~O1O1Means3000 && ~O2O1Means3000 && ~TAMeans3000 && ~PUCMeans3000  %E8 alone, means >3000 seconds
    if FitE8indiv || FitE8glbl
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, {R},{FitE8},{'E8'},{'E8'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, [],[],{'E8'});
    end
elseif E8Meanssub03000 && ~O1O1Meanssub03000 && ~O2O1Meanssub03000 && ~TAMeanssub03000 && ~PUCMeans3000 %E8 alone, meanssub0
    if FitE8indiv || FitE8glbl
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8sub0}, {SEsE8sub0},...
            {'Oid-E894-O1'}, {R},{FitE8},{'E8'},{'E8'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8sub0}, {SEsE8sub0},...
            {'Oid-E894-O1'}, [],[],{'E8'});
    end
elseif E8allmeans && TAallmeans  && ~O1O1allmeans && ~O2O1allmeans && ~PUCMeans3000  %E8 and TA, means
    if FitE8indiv || FitE8glbl
        [hE8,PlotHandleE8] = plotconccurve({concE8,concTA}, {pLoopsE8, pLoopsTA}, {SEsE8,SEsTA},...
            {'Oid-E894-O1','Oid-TA94-O1'}, {R,R},{FitE8,FitTA},{'E8','TA'},{'E8','TA'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, [],[],{'E8'});
    end
elseif E8Means3000 && TAMeans3000 && ~O1O1Means3000 && ~O2O1Means3000 && ~PUCMeans3000   %E8 and TA, means >3000 seconds
    if FitE8indiv || FitE8glbl
        %FIX FROM HERE
        [hE8,PlotHandleE8] = plotconccurve({concE8,concTA}, {pLoopsE8, pLoopsTA}, {SEsE8,SEsTA},...
            {'Oid-E894-O1','Oid-TA94-O1'}, {R,R},{FitE8,FitTA},{'E8','TA'},{'E8','TA'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8}, {pLoopsE8}, {SEsE8},...
            {'Oid-E894-O1'}, [],[],{'E8'});
    end
elseif (E8Meanssub03000 && TAMeanssub03000 && ~O1O1Meanssub03000 && ~O2O1Meanssub03000 && ~PUCMeans3000) || ...
        (E8Meanssub023000 && TAMeanssub023000 && ~O1O1Meanssub023000 && ~O2O1Meanssub023000 && ~PUCMeans3000)%E8 and TA, meanssub0 or sub02
    if FitTAindiv && FitTAglbl &&~FitE8glbl
        [hE8,PlotHandleE8] = plotconccurve({concE8,concTA}, {pLoopsE8sub0, pLoopsTAsub0}, {SEsE8sub0,SEsTAsub0},...
            {'Oid-E894-O1','Oid-TA94-O1'}, {R,R,R},{FitE8wTA,FitTAglobal,FitTA},{'E8','TA'},{'E8','TA','--r'});
    elseif FitE8glbl && FitTAindiv && FitTAglbl
        [hE8,PlotHandleE8] = plotconccurve({concE8,concTA}, {pLoopsE8sub0, pLoopsTAsub0}, {SEsE8sub0,SEsTAsub0},...
            {'Oid-E894-O1','Oid-TA94-O1'}, {R,R,R,R},{FitE8global,FitE8wTA,FitTAglobal,FitTA},{'E8','TA'},{':k','E8','TA','--r'});
    else
        [hE8,PlotHandleE8] = plotconccurve({concE8,concTA}, {pLoopsE8sub0, pLoopsTAsub0}, {SEsE8sub0,SEsTAsub0},...
            {'Oid-E894-O1','Oid-TA94-O1'}, [],[],{'E8','TA'});
    end
%Operator change plots
elseif E8Means3000 && O1O1Means3000 && O2O1Means3000 && ~TAMeans3000 && ~PUCMeans3000
    if FitE8indiv && ~FitE8glbl %assume if want to have one fit, want all fits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1}, {SEsE8,SEsE8O1O1,SEsE8O2O1},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, {R,R,R},{FitE8,FitO1O1,FitO2O1},...
            {'E8','E8O1O1','E8O2O1'},{'E8','E8O1O1','E8O2O1'});
    elseif FitE8glbl && ~FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1}, {SEsE8,SEsE8O1O1,SEsE8O2O1},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, {R,R,R},{FitE8global,FitO1O1global,FitO2O1global},...
            {'E8','E8O1O1','E8O2O1'},{'E8','E8O1O1','E8O2O1'});
    elseif FitE8glbl && FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1}, {SEsE8,SEsE8O1O1,SEsE8O2O1},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, ...
            {R,R,R,R,R,R},{FitE8, FitO1O1, FitO2O1,FitE8global,FitO1O1global,FitO2O1global},...
            {'E8','E8O1O1','E8O2O1'},{'--k','--b','--m','E8','E8O1O1','E8O2O1'});
    else %nofits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1}, {SEsE8,SEsE8O1O1,SEsE8O2O1},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, [],[],...
            {'E8','E8O1O1','E8O2O1'});
    end
%elseif (E8Meanssub03000 && O1O1Meanssub03000 && O2O1Meanssub03000 && ~TAMeanssub03000 && ~PUCMeans3000) || ...
        %(E8Meanssub023000 && O1O1Meanssub023000 && O2O1Meanssub023000 && ~TAMeanssub023000 && ~PUCMeans3000)
   %As of 3/2011 decided want to pair O2-O1 sub0 with sub02 for the rest
elseif (E8Meanssub03000 && O1O1Meanssub03000 && O2O1Meanssub03000 && ~TAMeanssub03000 && ~PUCMeans3000) || ...
        (E8Meanssub023000 && O1O1Meanssub023000 && O2O1Meanssub03000 && ~TAMeanssub023000 && ~PUCMeans3000)
    if FitE8indiv && ~FitE8glbl %assume if want to have one fit, want all fits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, {R,R,R},{FitE8,FitO1O1,FitO2O1},...
            {'E8','E8O1O1','E8O2O1'},{'E8','E8O1O1','E8O2O1'});
    elseif FitE8glbl && ~FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, {R,R,R},{FitE8global,FitO1O1global,FitO2O1global},...
            {'E8','E8O1O1','E8O2O1'},{'E8','E8O1O1','E8O2O1'});
    elseif FitE8glbl && FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, ...
            {R,R,R,R,R,R},{FitE8, FitO1O1, FitO2O1,FitE8global,FitO1O1global,FitO2O1global},...
            {'E8','E8O1O1','E8O2O1'},{'--k','--b','--m','E8','E8O1O1','E8O2O1'});
    else %nofits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1'}, [],[],...
            {'E8','E8O1O1','E8O2O1'});
    end
%Operator change, and TA
elseif E8Means3000 && O1O1Means3000 && O2O1Means3000 && TAMeans3000 && ~PUCMeans3000
    if FitE8indiv && ~FitE8glbl %assume if want to have one fit, want all fits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1,pLoopsTA}, {SEsE8,SEsE8O1O1,SEsE8O2O1,SEsTA},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, {R,R,R,R},{FitE8,FitO1O1,FitO2O1,FitTA},...
            {'E8','E8O1O1','E8O2O1','TA'},{'E8','E8O1O1','E8O2O1','TA'});
    elseif FitE8glbl && ~FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1,pLoopsTA}, {SEsE8,SEsE8O1O1,SEsE8O2O1,SEsTA},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, {R,R,R,R},{FitE8global,FitO1O1global,FitO2O1global,FitTAglobal},...
            {'E8','E8O1O1','E8O2O1','TA'},{'E8','E8O1O1','E8O2O1','TA'});
    elseif FitE8glbl && FitE8indiv %This plots three sets of curves: individual fits for all 4 data sets, plus the global
            %fit with just the E8 data sets, plus the global fit with all 4
            %data sets
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1,pLoopsTA}, {SEsE8,SEsE8O1O1,SEsE8O2O1,SEsTA},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, ...
            {R,R,R,R,R,R,R,R,R,R,R},{FitE8, FitO1O1, FitO2O1,FitTA, ...
            FitE8global,FitO1O1global,FitO2O1global, ...
            FitE8wTA,FitO1O1wTA,FitO2O1wTA,FitTAglobal},...
            {'E8','E8O1O1','E8O2O1','TA'},{':k',':b',':m',':r','--k','--b','--m','E8','E8O1O1','E8O2O1','TA'});
    else %nofits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8,pLoopsE8O1O1,pLoopsE8O2O1,pLoopsTA}, {SEsE8,SEsE8O1O1,SEsE8O2O1,SEsTA},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','TA'}, [],[],...
            {'E8','E8O1O1','E8O2O1','TA'});
    end
%elseif (E8Meanssub03000 && O1O1Meanssub03000 && O2O1Meanssub03000 && TAMeanssub03000 && ~PUCMeans3000) || ...
        %(E8Meanssub023000 && O1O1Meanssub023000 && O2O1Meanssub023000 && TAMeanssub023000 && ~PUCMeans3000) 
   %As of 3/2011 decided to use O2-O1 sub0 with sub02 for the rest
elseif (E8Meanssub03000 && O1O1Meanssub03000 && O2O1Meanssub03000 && TAMeanssub03000 && ~PUCMeans3000) || ...
        (E8Meanssub023000 && O1O1Meanssub023000 && O2O1Meanssub03000 && TAMeanssub023000 && ~PUCMeans3000) 
    if FitE8indiv && ~FitE8glbl %assume if want to have one fit, want all fits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0,pLoopsTAsub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, {R,R,R,R},{FitE8,FitO1O1,FitO2O1,FitTA},...
            {'E8','E8O1O1','E8O2O1','TA'},{'E8','E8O1O1','E8O2O1','TA'});
    elseif FitE8glbl && ~FitE8indiv
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0,pLoopsTAsub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, {R,R,R,R},{FitE8global,FitO1O1global,FitO2O1global,FitTAglobal},...
            {'E8','E8O1O1','E8O2O1','TA'},{'E8','E8O1O1','E8O2O1','TA'});
    elseif FitE8glbl && FitE8indiv %This plots three sets of curves: individual fits for all 4 data sets, plus the global
            %fit with just the E8 data sets, plus the global fit with all 4
            %data sets
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0,pLoopsTAsub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'}, ...
            {R,R,R,R,R,R,R,R,R,R,R},{FitE8, FitO1O1, FitO2O1,FitTA, ...
            FitE8global,FitO1O1global,FitO2O1global, ...
            FitE8wTA,FitO1O1wTA,FitO2O1wTA,FitTAglobal},...
            {'E8','E8O1O1','E8O2O1','TA'},{':k',':b',':m',':r','--k','--b','--m','E8','E8O1O1','E8O2O1','TA'});
    else %nofits
        [hE8O1O2,PlotHandleE8O1O2] = plotconccurve({concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsE8sub0,pLoopsE8O1O1sub0,pLoopsE8O2O1sub0,pLoopsTAsub0}, {SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'Oid-E894-O1','O1-E894-O1','O2-E894-O1','TA'}, [],[],...
            {'E8','E8O1O1','E8O2O1','TA'});
    end
elseif PUCMeans3000 && PUCBMeans3000 && PUCMMeans3000 && E8Means3000 && O1O1Means3000 && ...
        O2O1Means3000 && TAMeans3000 %Oid-/O1-/O2-E8, and TA, and PUC
    if FitE8TAPUC
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8, pLoopsE8O1O1, pLoopsE8O2O1, pLoopsTA},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8,SEsE8O1O1,SEsE8O2O1,SEsTA},...
            {'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'},...
            {R,R,R,R,R,R,R},{FitPUCtotwET,FitPUCMwET,FitPUCBwET,FitE8wPUC,FitO1O1wPUC,FitO2O1wPUC,FitTAwPUC},...
            {'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'},{'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'});
    else %no fits 
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8sub0, pLoopsE8O1O1sub0, pLoopsE8O2O1sub0, pLoopsTAsub0},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'},...
            [],[],{'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'});
    end    
elseif (PUCMeans3000 && PUCBMeans3000 && PUCMMeans3000 && E8Meanssub03000 && O1O1Meanssub03000 && ...
        O2O1Meanssub03000 && TAMeanssub03000) || (PUCMeans3000 && PUCBMeans3000 && PUCMMeans3000 && ...
        E8Meanssub023000 && O1O1Meanssub023000 && O2O1Meanssub03000 && TAMeanssub023000) %Oid-/O1-/O2-E8, and TA, sub0 or sub02, and PUC
    if FitE8TAPUCsub0 || FitE8TAPUCsub02
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8sub0, pLoopsE8O1O1sub0, pLoopsE8O2O1sub0, pLoopsTAsub0},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'},...
            {R,R,R,R,R,R,R},{FitPUCtotwET,FitPUCMwET,FitPUCBwET,FitE8wPUC,FitO1O1wPUC,FitO2O1wPUC,FitTAwPUC},...
            {'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'},{'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'});
    else %no fits 
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8,concE8O1O1,concE8O2O1,concTA}, ...
            {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8sub0, pLoopsE8O1O1sub0, pLoopsE8O2O1sub0, pLoopsTAsub0},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8sub0,SEsE8O1O1sub0,SEsE8O2O1sub0,SEsTAsub0},...
            {'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1','O1-E894-O1','O2-E894-O1','Oid-TA94-O1'},...
            [],[],{'PUC','PUCM','PUCB','E8','E8O1O1','E8O2O1','TA'});
    end
elseif PUCMeans3000 && ~PUCBMeans3000 && E8Meanssub03000 %E8 and PUC total only
    if FitPUCtotal
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concE8}, {pLoopsPUC,pLoopsE8sub0},...
            {SEsPUC,SEsE8sub0},{'O1-PUC306-Oid','Oid-E894-O1'},...
            {R,R},{FitPUC,FitE8global},{'PUC','E8'},{'PUC','E8'});
    else
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concE8}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R},{FitPUC,FitE8global},{'PUC','E8'},{'PUC','E8'});
    end
        
elseif PUCMeans3000 && PUCBMeans3000 && PUCMMeans3000 && (E8Meanssub03000 || E8Meanssub023000) %E8 and PUC with the looped states separated out, beads more than 3000 seconds
    if FitPUCtotal
        %FIX
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concE8}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R},{FitPUC,FitE8global},{'PUC','E8'},{'PUC','E8'});
    elseif FitPUCBvsM || FitPUCBvsMsameK
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8sub0},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8sub0},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1'},...
            {R,R,R,R},{FitPUCtot,FitPUCM,FitPUCB,FitE8global},{'PUC','PUCM','PUCB','E8'},{'PUC','PUCM','PUCB','E8'});
    else %no fits 
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC, concE8}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB, pLoopsE8sub0},...
            {SEsPUC,SEsPUCM,SEsPUCB,SEsE8sub0},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B','Oid-E894-O1'},...
            [],[],{'PUC','PUCM','PUCB','E8'});
    end
elseif PUCallmeans && PUCBmeans && PUCMmeans %PUC with the looped states separated out
    if FitPUCtotal
        %FIX
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R,R},{FitPUC,FitPUCM,FitPUCB},{'PUC','PUCM','PUCB'},{'PUC','PUCM','PUCB'});
    elseif FitPUCBvsM
        %FIX
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R,R},{FitPUC,FitPUCM,FitPUCB},{'PUC','PUCM','PUCB'},{'PUC','PUCM','PUCB'});
    else %no fits 
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            [],[],{'PUC','PUCM','PUCB'});
    end
elseif PUCMeans3000 && PUCBMeans3000 && PUCMMeans3000 %PUC with the looped states separated out, beads more than 3000 seconds
    if FitPUCtotal
        %FIX
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R,R},{FitPUCtot,FitPUCM,FitPUCB},{'PUC','PUCM','PUCB'},{'PUC','PUCM','PUCB'});
    elseif FitPUCBvsM || FitPUCBvsMsameK
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            {R,R,R},{FitPUCtot,FitPUCM,FitPUCB},{'PUC','PUCM','PUCB'},{'PUC','PUCM','PUCB'});
    else %no fits 
        [hPUC,PlotHandlePUC] = plotconccurve({concPUC,concPUC,concPUC}, {pLoopsPUC,pLoopsPUCM,pLoopsPUCB},...
            {SEsPUC,SEsPUCM,SEsPUCB},{'O1-PUC306-Oid, B+M','O1-PUC306-Oid, M','O1-PUC306-Oid, B'},...
            [],[],{'PUC','PUCM','PUCB'});
    end
    
end
%% FITS: perform fits to obtain fit parameters
R=logspace(-15,-5,1000);
%% E8 indiv fit
[concE8, pLoopsDistribsE8] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',3000);
%[KidE8,K1E8,JE8,SEKidE8,SEK1E8,SEJE8] = fitconccurve(concE8.*10^12,pLoopsDistribsE8);
[Kid,K1,JE8,SEKid,SEK1,SEJ] = fitconccurve2(concE8.*10^12,pLoopsDistribsE8,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfit.mat');
%FitE8 = plotpLoopTheory(KidE8*10^-12,K1E8*10^-12,JE8*10^-12);
%close
clear concE8 pLoopsDistribsE8 Kid K1 J SEKid SEK1 SEJ

%% E8 indiv fit sub0
[concE8, pLoopsDistribsE8sub0] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
%[KidE8sub0,K1E8sub0,JE8sub0,SEKidE8sub0,SEK1E8sub0,SEJE8sub0] = fitconccurve(concE8.*10^12,pLoopsDistribsE8sub0);
[Kid,K1,JE8,SEKid,SEK1,SEJ] = fitconccurve2(concE8.*10^12,pLoopsDistribsE8sub0,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfitsub0.mat');
%FitE8 = plotpLoopTheory(KidE8sub0*10^-12,K1E8sub0*10^-12,JE8sub0*10^-12);
%close
clear concE8 pLoopsDistribsE8sub0 KidE8sub0 K1E8sub0 JE8sub0 SEKidE8sub0 SEK1E8sub0 SEJE8sub0

%% E8 indiv fit sub02
[concE8, pLoopsDistribsE8sub02] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
%[KidE8sub0,K1E8sub0,JE8sub0,SEKidE8sub0,SEK1E8sub0,SEJE8sub0] = fitconccurve(concE8.*10^12,pLoopsDistribsE8sub0);
[Kid,K1,JE8,SEKid,SEK1,SEJ] = fitconccurve2(concE8.*10^12,pLoopsDistribsE8sub02,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/E894 conc curve/110216E8indivfitsub02.mat');
%FitE8 = plotpLoopTheory(KidE8sub0*10^-12,K1E8sub0*10^-12,JE8sub0*10^-12);
%close
clear concE8 pLoopsDistribsE8sub0 KidE8sub0 K1E8sub0 JE8sub0 SEKidE8sub0 SEK1E8sub0 SEJE8sub0

%% TA indiv fit
[concTA,pLoopsDistribsTA] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
%[KidTA,K1TA,JTA,SEKidTA,SEK1TA,SEJTA] = fitconccurve(concTA.*10^12,pLoopsDistribsTA);
[Kid,K1,JTA,SEKid,SEK1,SEJTA] = fitconccurve2(concTA.*10^12,pLoopsDistribsTA,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfit.mat');
%FitTA = plotpLoopTheory(KidTA,K1TA,JTA);
%close
%save('/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/TAindivfitparams','KidTA','K1TA','JTA','SEKidTA','SEK1TA','SEJTA');
clear concTA pLoopsDistribsTA Kid K1 JTA SEKid SEK1 SEJTA

%% TA indiv fit, sub0
[concTA,pLoopsDistribsTAsub0] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
%[K1TAsub0,K2TAsub0,JTAsub0,SEK1TAsub0,SEK2TAsub0,SEJTAsub0] = fitconccurve(concTA.*10^12,pLoopsDistribsTAsub0);
[Kid,K1,JTA,SEKid,SEK1,SEJTA] = fitconccurve2(concTA.*10^12,pLoopsDistribsTAsub0,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfitsub0.mat');
%FitTAsub0 = plotpLoopTheory(KidTAsub0,K1TAsub0,JTAsub0);
%close
clear concTA pLoopsDistribsTAsub0 K1TAsub0 K2TAsub0 JTAsub0 SEK1TAsub0 SEK2TAsub0 SEJTAsub0
%FitTAindivsub0 = 1;

%% TA indiv fit, sub02
[concTA,pLoopsDistribsTAsub02] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
%[K1TAsub0,K2TAsub0,JTAsub0,SEK1TAsub0,SEK2TAsub0,SEJTAsub0] = fitconccurve(concTA.*10^12,pLoopsDistribsTAsub0);
[Kid,K1,JTA,SEKid,SEK1,SEJTA] = fitconccurve2(concTA.*10^12,pLoopsDistribsTAsub02,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/TA94 conc curve/110216TAindivfitsub02.mat');
%FitTAsub0 = plotpLoopTheory(KidTAsub0,K1TAsub0,JTAsub0);
%close
clear concTA pLoopsDistribsTAsub02 K1TAsub0 K2TAsub0 JTAsub0 SEK1TAsub0 SEK2TAsub0 SEJTAsub0
%FitTAindivsub0 = 1;

%% O1E894O1 indiv fit
[concE8O1O1,dataE8O1O1] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
%[K1O1O1,K1O1O1,JO1O1,SEK1O1O1,SEK1O1O1,SEJO1O1] = fitconccurve(concE8O1O1.*10^12,dataE8O1O1,'samek');
[K1,K1,J,SEK1,SEK1,SEJ] = fitconccurve2(concE8O1O1.*10^12,dataE8O1O1,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfit.mat','sameK');
%FitE8O1O1 = plotpLoopTheory(K1O1O1,K1O1O1,JO1O1);
%close
clear concE8O1O1 dataE8O1O1 K1 K1 J SEK1 SEK1 SEJ
%FitO1O1indiv = 1;

%% O1E894O1 indiv fit, sub0
[concE8O1O1,dataE8O1O1sub0] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
%[K1O1O1sub0,K1O1O1sub0,JO1O1sub0,SEK1O1O1sub0,SEK1O1O1sub0,SEJO1O1sub0] = fitconccurve(concE8O1O1.*10^12,dataE8O1O1sub0,'samek');
[K1,K1,J,SEK1,SEK1,SEJ] = fitconccurve2(concE8O1O1.*10^12,dataE8O1O1sub0,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfitsub0.mat','sameK');
%FitE8O1O1sub0 = plotpLoopTheory(K1O1O1sub0,K1O1O1sub0,JO1O1sub0);
%close
clear concE8O1O1 dataE8O1O1sub0 K1 K1 J SEK1 SEK1 SEJ
%FitO1O1sub0indiv = 1;

%% O1E894O1 indiv fit, sub02
[concE8O1O1,dataE8O1O1sub02] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
%[K1O1O1sub0,K1O1O1sub0,JO1O1sub0,SEK1O1O1sub0,SEK1O1O1sub0,SEJO1O1sub0] = fitconccurve(concE8O1O1.*10^12,dataE8O1O1sub0,'samek');
[K1,K1,J,SEK1,SEK1,SEJ] = fitconccurve2(concE8O1O1.*10^12,dataE8O1O1sub02,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O1E894O1 conc curve/110216O1O1indivfitsub02.mat','sameK');
%FitE8O1O1sub0 = plotpLoopTheory(K1O1O1sub0,K1O1O1sub0,JO1O1sub0);
%close
clear concE8O1O1 dataE8O1O1sub0 K1 K1 J SEK1 SEK1 SEJ

%% O2E894O1 indiv fit
[concE8O2O1,dataE8O2O1] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
%[K1O2O1,K2O2O1,JO2O1,SEK1O2O1,SEK2O2O1,SEJO2O1] = fitconccurve(concE8O2O1.*10^12,dataE8O2O1);
[K1,K2,J,SEK1,SEK2,SEJ] = fitconccurve2(concE8O2O1.*10^12,dataE8O2O1,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfit.mat');
%FitE8O1O1 = plotpLoopTheory(K1O1O1,K1O1O1,JO1O1);
%close
clear concE8O1O1 dataE8O1O1 K1 K2 J SEK1 SEK2 SEJ
%FitO1O1indiv = 1;

%% O2E894O1 indiv fit, sub0
[concE8O2O1,dataE8O2O1sub0] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
%[K1O2O1sub0,K2O2O1sub0,JO2O1sub0,SEK1O2O1sub0,SEK2O2O1sub0,SEJO2O1sub0] = fitconccurve(concE8O2O1.*10^12,dataE8O2O1sub0);
[K1,K2,J,SEK1,SEK2,SEJ] = fitconccurve2(concE8O2O1.*10^12,dataE8O2O1sub0,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfitsub0.mat');
%FitE8O2O1sub0 = plotpLoopTheory(K1O2O1sub0,K2O2O1sub0,JO2O1sub0);
%close
clear concE8O2O1 dataE8O2O1sub0 K1 K2 J SEK1 SEK2 SEJ
%FitO2O1sub0indiv = 1;

%% O2E894O1 indiv fit, sub02
[concE8O2O1,dataE8O2O1sub02] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
%[K1O2O1sub0,K2O2O1sub0,JO2O1sub0,SEK1O2O1sub0,SEK2O2O1sub0,SEJO2O1sub0] = fitconccurve(concE8O2O1.*10^12,dataE8O2O1sub0);
[K1,K2,J,SEK1,SEK2,SEJ] = fitconccurve2(concE8O2O1.*10^12,dataE8O2O1sub02,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/O2E894O1 conc curve/110216O2O1indivfitsub02.mat');
%FitE8O2O1sub0 = plotpLoopTheory(K1O2O1sub0,K2O2O1sub0,JO2O1sub0);
%close
clear concE8O2O1 dataE8O2O1sub02 K1 K2 J SEK1 SEK2 SEJ
%FitO2O1sub0indiv = 1;

%% E8OidO1O2 global fits
[fitparams, paramSEs] = fitconccurve_global('E8OidO1O2',0);
JE8 = fitparams(1); Kid = fitparams(2); K1 = fitparams(3); K2=fitparams(4);
SEJE8 = paramSEs(1); SEKid = paramSEs(2); SEK1 = paramSEs(3); SEK2 = paramSEs(4);

%% E8OidO1O2sub0 global fit
[fitparams, paramSEs] = fitconccurve_global('E8OidO1O2sub0',0);
JE8 = fitparams(1); Kid = fitparams(2); K1 = fitparams(3); K2=fitparams(4);
SEJE8 = paramSEs(1); SEKid = paramSEs(2); SEK1 = paramSEs(3); SEK2 = paramSEs(4);

%% PUC total (nonglobal) fit
[concPUC, pLoopsDistribsPUC] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribs',3000);
[K1,K2,J,SEK1,SEK2,SEJ] = fitconccurve2(concPUC.*10^12,pLoopsDistribsPUC,...
    '/Volumes/dumbo-3/stephj/TPM data analysis/SJLacI/PUC306 conc curve/110218PUCtotfit.mat');
%FitPUC = plotpLoopTheory(KidPUC*10^-12,K1PUC*10^-12,JPUC*10^-12);
%close
clear concPUC pLoopsDistribsPUC KidPUC K1PUC JPUC SEKidPUC SEK1PUC SEJPUC

%% PUC B&M global fit
[fitparams, paramSEs] = fitconccurve_global('PUCBvsM',0);
JB = fitparams(1); JM = fitparams(2); Kid = fitparams(3); K1 = fitparams(4);
SEJB = paramSEs(1); SEJM = paramSEs(2); SEKid = paramSEs(3); SEK1 = paramSEs(4);

%% PUC B&M global fit, enforce E8 Kd's (Kid = 9.3 pM, K1 = 42.9 pM)
[fitparams,paramSEs]=fitconccurve_global('PUCBvsM',0); %Hand-coded the Kd's in temporarily
%Result was JB = 392.9738 +/- 37.1741, JM = 478.1499 +/- 45.8122

%% E107BvsM fit with Kd's from E8 and TA conc curves:
[fitparams,paramSEs]=fitconccurve_global('E8107BvsMsub02',0);



