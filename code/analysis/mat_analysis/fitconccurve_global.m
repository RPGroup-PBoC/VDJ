%[fitparams,paramSEs]=fitconccurve_global(dataset,dimers)
%
%For doing a global fit to multiple concentration curves.  Based off of
%globalfit.m.
%
%Inputs: "dimers" is either 1 or 0--IMPLEMENTED FOR E8OidO1O2TA(/sub0/sub02) and E8OidO1O2TAPUC(/sub0/sub02) only!  "Dataset" can be
%'E8OidO1O2','E8OidO1O2sub0','E8OidO1O2TA', 'E8OidO1O2TAsub0' or 'PUCBvsM'.
%
%10/2011 added a way to fit a new concentration curve with the Kd's from
%the E8OidO1O2TAsub02 global fit, so that I can fit E8107BvsM (to use this
%option dataset = 'E8107BvsMsub02').
%
%As a reminder: "model" is of the form @(params)[(eqn1-data1)./SEs1; (eqn1-data1)./SEs1; etc]
%
%Steph 12/10

function [fitparams,paramSEs]=fitconccurve_global(dataset,dimers)

%As with fitconccurve and fitconccurve/_general, the errors will be
%bootstrapped:
nboot = 10000;
nbins = nboot/100;

%Set fit options
opt=optimset(@lsqnonlin);
opt=optimset(opt,'tolx',1e-12,'tolfun',1e-12,'maxfunevals',1e5,'maxiter',1e5);

if strcmpi(dataset,'E8OidO1O2')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2sub0')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2sub02')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %[concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %As of 3/2011 decided that the method used in sub0, not sub02, was the
    %best for O2-O1 (see pLoopsHistosData tex files for details about the
    %differences between the two treatements)
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2TA')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2TAsub0')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;

elseif strcmpi(dataset,'E8OidO1O2TAsub02')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %[concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %As of 3/2011 decided that the method used in sub0, not sub02, was the
    %best for O2-O1 (see pLoopsHistosData tex files for details about the
    %differences between the two treatements)
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;
    
elseif strcmpi(dataset,'PUCBvsM')
    [concPUC, pLoopsPUCBdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsB',3000);
    [concPUC, pLoopsPUCMdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsM',3000);
    
    [newmeansPUCB,newSEsPUCB,pLoopsPUCB,SEsPUCB] = bootstrap_concs(concPUC,pLoopsPUCBdistrib);
    [newmeansPUCM,newSEsPUCM,pLoopsPUCM,SEsPUCM] = bootstrap_concs(concPUC,pLoopsPUCMdistrib);

    concPUC = concPUC.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2TAPUC')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribs',3000);
    [concPUC, pLoopsPUCBdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsB',3000);
    [concPUC, pLoopsPUCMdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsM',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);
    [newmeansPUCB,newSEsPUCB,pLoopsPUCB,SEsPUCB] = bootstrap_concs(concPUC,pLoopsPUCBdistrib);
    [newmeansPUCM,newSEsPUCM,pLoopsPUCM,SEsPUCM] = bootstrap_concs(concPUC,pLoopsPUCMdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;
    concPUC = concPUC.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2TAPUCsub0')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concPUC, pLoopsPUCBdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsB',3000);
    [concPUC, pLoopsPUCMdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsM',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);
    [newmeansPUCB,newSEsPUCB,pLoopsPUCB,SEsPUCB] = bootstrap_concs(concPUC,pLoopsPUCBdistrib);
    [newmeansPUCM,newSEsPUCM,pLoopsPUCM,SEsPUCM] = bootstrap_concs(concPUC,pLoopsPUCMdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;
    concPUC = concPUC.*10^12;
    
elseif strcmpi(dataset,'E8OidO1O2TAPUCsub02')
    [concOidO1, pLoopsOidO1distrib] = LoadDataAnalysis('E894 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    [concO1O1, pLoopsO1O1distrib] = LoadDataAnalysis('O1E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %[concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    %As of 3/2011 decided that the method used in sub0, not sub02, was the
    %best for O2-O1 (see pLoopsHistosData tex files for details about the
    %differences between the two treatements)
    [concO2O1, pLoopsO2O1distrib] = LoadDataAnalysis('O2E894O1 conc curve','HistoAnal','SJLacI','alldistribssub0',3000);
    [concTA, pLoopsTAdistrib] = LoadDataAnalysis('TA94 conc curve','HistoAnal','SJLacI','alldistribssub02',3000);
    [concPUC, pLoopsPUCBdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsB',3000);
    [concPUC, pLoopsPUCMdistrib] = LoadDataAnalysis('PUC306 conc curve','HistoAnal','SJLacI','alldistribsM',3000);
    
    [newmeansOidO1,newSEsOidO1,pLoopsOidO1,SEsOidO1] = bootstrap_concs(concOidO1,pLoopsOidO1distrib);
    [newmeansO1O1,newSEsO1O1,pLoopsO1O1,SEsO1O1] = bootstrap_concs(concO1O1,pLoopsO1O1distrib);
    [newmeansO2O1,newSEsO2O1,pLoopsO2O1,SEsO2O1] = bootstrap_concs(concO2O1,pLoopsO2O1distrib);
    [newmeansTA,newSEsTA,pLoopsTA,SEsTA] = bootstrap_concs(concTA,pLoopsTAdistrib);
    [newmeansPUCB,newSEsPUCB,pLoopsPUCB,SEsPUCB] = bootstrap_concs(concPUC,pLoopsPUCBdistrib);
    [newmeansPUCM,newSEsPUCM,pLoopsPUCM,SEsPUCM] = bootstrap_concs(concPUC,pLoopsPUCMdistrib);

    concOidO1 = concOidO1.*10^12; %Since I do all my fitting in pM units
    concO1O1 = concO1O1.*10^12;
    concO2O1 = concO2O1.*10^12;
    concTA = concTA.*10^12;
    concPUC = concPUC.*10^12;
    
elseif strcmpi(dataset,'E8107BvsMsub02')
    [concE107, pLoopsE107distrib] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','alldistribssub02',3000);
    [concE107M, pLoopsE107distribM] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','alldistribssub02M',3000);
    [concE107B, pLoopsE107distribB] = LoadDataAnalysis('E8107conccurve','ThreshAnal','SJLacI','alldistribssub02B',3000);
    
    [newmeansE107,newSEsE107,pLoopsE107,SEsE107] = bootstrap_concs(concE107,pLoopsE107distrib);
    [newmeansE107M,newSEsE107M,pLoopsE107M,SEsE107M] = bootstrap_concs(concE107M,pLoopsE107distribM);
    [newmeansE107B,newSEsE107B,pLoopsE107B,SEsE107B] = bootstrap_concs(concE107B,pLoopsE107distribB);

    concE107= concE107.*10^12; %Since I do all my fitting in pM units
    concE107M = concE107M.*10^12;
    concE107B = concE107B.*10^12;
    
    %Load the Kd's to enforce:
    temp = load('/Volumes/dumbo-4/stephj/TPM data analysis/SJLacI/110309E8OidO1O2TAfitsub02.mat');

    Kid = temp.fitparams(3);
    K1 = temp.fitparams(4);
    Kiddistrib = temp.paramsbs(:,3);
    K1distrib = temp.paramsbs(:,4);
    
end

missed = 0;

if strcmpi(dataset,'E8OidO1O2') || strcmpi(dataset,'E8OidO1O2sub0') || strcmpi(dataset,'E8OidO1O2sub02')
    %params is [J, Kid, K1, K2]
    paramnames{1} = 'J';
    paramnames{2} = 'Kid';
    paramnames{3} = 'K1';
    paramnames{4} = 'K2';

    %Since I want to be able to have different numbers of
    %concentrations for each data set, I need to make a separate
    %function for each data set:  (is that really true?)

    f_OidO1 = @(params)(((0.5.*concOidO1.*params(1)./(params(2)*params(3)))./(1+...
        concOidO1./params(2)+concOidO1./params(3)+concOidO1.^2./(params(2)*params(3))+...
        0.5.*concOidO1.*params(1)./(params(2)*params(3)))-pLoopsOidO1)./SEsOidO1);
    f_O1O1 = @(params)(((0.5.*concO1O1.*params(1)./(params(3)^2))./(1+...
        2*concO1O1./params(3)+concO1O1.^2./(params(3)^2)+...
        0.5.*concO1O1.*params(1)./(params(3)^2))-pLoopsO1O1)./SEsO1O1);
    f_O2O1 = @(params)(((0.5.*concO2O1.*params(1)./(params(3)*params(4)))./(1+...
        concO2O1./params(3)+concO2O1./params(4)+concO2O1.^2./(params(3)*params(4))+...
        0.5.*concO2O1.*params(1)./(params(3)*params(4)))-pLoopsO2O1)./SEsO2O1);

    model=@(params)([f_OidO1(params) f_O1O1(params) f_O2O1(params)]);


    %Use lsqnonlin to fit.  Second input are starting parameters, output is
    %best fit results for the params vector.  I can also add additional arguments
    %lower bound (vector), upper bound (vector) for all parameters.

    fitparams = lsqnonlin(model,[300,6,40,200],[0 0 0 0],[],opt);


    %Find the errors:
    paramsbs = zeros(nboot,length(fitparams));
    for j=1:nboot

        %If any of the SEs are 0, that's a problem for the fitting
        %routine (div by 0!)  O2's SEs can be 0 for the lowest
        %concentration.  So, same as in fitconccurve: make the SE for
        %any zero points to be the average of the rest of the SEs:
        if length(find(newSEsOidO1(j,:)))~=length(newSEsOidO1(j,:))
            newSEsOidO1(j,~logical(newSEsOidO1(j,:)))=mean(newSEsOidO1(j,logical(newSEsOidO1(j,:))));
        end
        if length(find(newSEsO1O1(j,:)))~=length(newSEsO1O1(j,:))
            newSEsO1O1(j,~logical(newSEsO1O1(j,:)))=mean(newSEsO1O1(j,logical(newSEsO1O1(j,:))));
        end
        if length(find(newSEsO2O1(j,:)))~=length(newSEsO2O1(j,:))
            newSEsO2O1(j,~logical(newSEsO2O1(j,:)))=mean(newSEsO2O1(j,logical(newSEsO2O1(j,:))));
        end



        newf_OidO1 = @(newparams)(((0.5.*concOidO1.*newparams(1)./(newparams(2)*newparams(3)))./(1+...
            concOidO1./newparams(2)+concOidO1./newparams(3)+concOidO1.^2./(newparams(2)*newparams(3))+...
            0.5.*concOidO1.*newparams(1)./(newparams(2)*newparams(3)))-newmeansOidO1(j,:))./newSEsOidO1(j,:));
        newf_O1O1 = @(newparams)(((0.5.*concO1O1.*newparams(1)./(newparams(3)^2))./(1+...
            2*concO1O1./newparams(3)+concO1O1.^2./(newparams(3)^2)+...
            0.5.*concO1O1.*newparams(1)./(newparams(3)^2))-newmeansO1O1(j,:))./newSEsO1O1(j,:));
        newf_O2O1 = @(newparams)(((0.5.*concO2O1.*newparams(1)./(newparams(3)*newparams(4)))./(1+...
            concO2O1./newparams(3)+concO2O1./newparams(4)+concO2O1.^2./(newparams(3)*newparams(4))+...
            0.5.*concO2O1.*newparams(1)./(newparams(3)*newparams(4)))-newmeansO2O1(j,:))./newSEsO2O1(j,:));

        newmodel=@(newparams)([newf_OidO1(newparams) newf_O1O1(newparams) newf_O2O1(newparams)]);

        try
            paramsbs(j,:) = lsqnonlin(newmodel,[300,6,40,200],[0 0 0 0],[],opt);
        catch
            missed = missed+1;
            oldparamsbs = paramsbs(1:j-1,:);
            clear paramsbs
            paramsbs=zeros(nboot-1,length(fitparams)); %I can't make the entries NaN's if the fitting fails because
                %std will be NaN then; so just remove that whole
                %attempt
            paramsbs(1:j-1,:) = oldparamsbs;
            keyboard
        end

        clear newmodel

%             if rem(j,10)==0
%                 disp(j)
%             end

    end
    paramSEs = std(paramsbs,1,1);
        
elseif strcmpi(dataset,'E8OidO1O2TA') || strcmpi(dataset,'E8OidO1O2TAsub0') || strcmpi(dataset,'E8OidO1O2TAsub02')
    %params is [JE8, JTA, Kid, K1,K2]
    paramnames{1} = 'JE8';
    paramnames{2} = 'JTA';
    paramnames{3} = 'Kid';
    paramnames{4} = 'K1';
    paramnames{5} = 'K2';

    %Since I want to be able to have different numbers of
    %concentrations for each data set, I need to make a separate
    %function for each data set:  (is that really true?)

    if dimers
        paramnames{6} = 'KDT';
        
        f_OidO1 = @(params)(((0.5.*concOidO1.*params(1).*(1 + ...
            params(6)./(8*concOidO1) - (1/8)*sqrt((params(6)./concOidO1).^2 + ...
            16*params(6)./concOidO1))./(params(3)*params(4)))./(1+...
            concOidO1./params(3)+concOidO1./params(4)+concOidO1.^2./(params(3)*params(4))+...
            0.5.*concOidO1.*params(1).*(1 + params(6)./(8*concOidO1) - ...
            (1/8)*sqrt((params(6)./concOidO1).^2 + 16*params(6)./concOidO1))./(params(3)*params(4))) - ...
            pLoopsOidO1)./SEsOidO1);
        f_O1O1 = @(params)(((0.5.*concO1O1.*params(1).*(1 + ...
            params(6)./(8*concO1O1) - (1/8)*sqrt((params(6)./concO1O1).^2 + ...
            16*params(6)./concO1O1))./(params(4)^2))./(1+...
            2*concO1O1./params(4)+concO1O1.^2./(params(4)^2)+...
            0.5.*concO1O1.*params(1).*(1 + params(6)./(8*concO1O1) - ...
            (1/8)*sqrt((params(6)./concO1O1).^2 + 16*params(6)./concO1O1))./(params(4)^2)) - ...
            pLoopsO1O1)./SEsO1O1);
        f_O2O1 = @(params)(((0.5.*concO2O1.*params(1).*(1 + ...
            params(6)./(8*concO2O1) - (1/8)*sqrt((params(6)./concO2O1).^2 + ...
            16*params(6)./concO2O1))./(params(4)*params(5)))./(1+...
            concO2O1./params(4)+concO2O1./params(5)+concO2O1.^2./(params(4)*params(5))+...
            0.5.*concO2O1.*params(1).*(1 + params(6)./(8*concO2O1) - ...
            (1/8)*sqrt((params(6)./concO2O1).^2 + 16*params(6)./concO2O1))./(params(4)*params(5))) - ...
            pLoopsO2O1)./SEsO2O1);
        f_TA = @(params)(((0.5.*concTA.*params(2).*(1 + ...
            params(6)./(8*concTA) - (1/8)*sqrt((params(6)./concTA).^2 + ...
            16*params(6)./concTA))./(params(3)*params(4)))./(1+...
            concTA./params(3)+concTA./params(4)+concTA.^2./(params(3)*params(4))+...
            0.5.*concTA.*params(2).*(1 + params(6)./(8*concTA) - ...
            (1/8)*sqrt((params(6)./concTA).^2 + 16*params(6)./concTA))./(params(3)*params(4))) - ...
            pLoopsTA)./SEsTA);
    else
        f_OidO1 = @(params)(((0.5.*concOidO1.*params(1)./(params(3)*params(4)))./(1+...
            concOidO1./params(3)+concOidO1./params(4)+concOidO1.^2./(params(3)*params(4))+...
            0.5.*concOidO1.*params(1)./(params(3)*params(4)))-pLoopsOidO1)./SEsOidO1);
        f_O1O1 = @(params)(((0.5.*concO1O1.*params(1)./(params(4)^2))./(1+...
            2*concO1O1./params(4)+concO1O1.^2./(params(4)^2)+...
            0.5.*concO1O1.*params(1)./(params(4)^2))-pLoopsO1O1)./SEsO1O1);
        f_O2O1 = @(params)(((0.5.*concO2O1.*params(1)./(params(4)*params(5)))./(1+...
            concO2O1./params(4)+concO2O1./params(5)+concO2O1.^2./(params(4)*params(5))+...
            0.5.*concO2O1.*params(1)./(params(4)*params(5)))-pLoopsO2O1)./SEsO2O1);
        f_TA = @(params)(((0.5.*concTA.*params(2)./(params(3)*params(4)))./(1+...
            concTA./params(3)+concTA./params(4)+concTA.^2./(params(3)*params(4))+...
            0.5.*concTA.*params(2)./(params(3)*params(4)))-pLoopsTA)./SEsTA);
    end

    model=@(params)([f_OidO1(params) f_O1O1(params) f_O2O1(params) f_TA(params)]);


    %Use lsqnonlin to fit.  Second input are starting parameters, output is
    %best fit results for the params vector.  I can also add additional arguments
    %lower bound (vector), upper bound (vector) for all parameters.

    if dimers
        [fitparams] = lsqnonlin(model,[300,1000,6,40,200,0.1],[0 0 0 0 0 0],[],opt);
    else
        [fitparams] = lsqnonlin(model,[300,1000,6,40,200],[0 0 0 0 0],[],opt);
    end


    %Find the errors:
    paramsbs = zeros(nboot,length(fitparams));
    for j=1:nboot

        %If any of the SEs are 0, that's a problem for the fitting
        %routine (div by 0!)  O2's SEs can be 0 for the lowest
        %concentration.  So, same as in fitconccurve: make the SE for
        %any zero points to be the average of the rest of the SEs:
        if length(find(newSEsOidO1(j,:)))~=length(newSEsOidO1(j,:))
            newSEsOidO1(j,~logical(newSEsOidO1(j,:)))=mean(newSEsOidO1(j,logical(newSEsOidO1(j,:))));
        end
        if length(find(newSEsO1O1(j,:)))~=length(newSEsO1O1(j,:))
            newSEsO1O1(j,~logical(newSEsO1O1(j,:)))=mean(newSEsO1O1(j,logical(newSEsO1O1(j,:))));
        end
        if length(find(newSEsO2O1(j,:)))~=length(newSEsO2O1(j,:))
            newSEsO2O1(j,~logical(newSEsO2O1(j,:)))=mean(newSEsO2O1(j,logical(newSEsO2O1(j,:))));
        end
        if length(find(newSEsTA(j,:)))~=length(newSEsTA(j,:))
            newSEsTA(j,~logical(newSEsTA(j,:)))=mean(newSEsTA(j,logical(newSEsTA(j,:))));
        end

        if dimers
            newf_OidO1 = @(newparams)(((0.5.*concOidO1.*newparams(1).*(1 + ...
                newparams(6)./(8*concOidO1) - (1/8)*sqrt((newparams(6)./concOidO1).^2 + ...
                16*newparams(6)./concOidO1))./(newparams(3)*newparams(4)))./(1+...
                concOidO1./newparams(3)+concOidO1./newparams(4)+concOidO1.^2./(newparams(3)*newparams(4))+...
                0.5.*concOidO1.*newparams(1).*(1 + ...
                newparams(6)./(8*concOidO1) - (1/8)*sqrt((newparams(6)./concOidO1).^2 + ...
                16*newparams(6)./concOidO1))./(newparams(3)*newparams(4)))-newmeansOidO1(j,:))./newSEsOidO1(j,:));
            newf_O1O1 = @(newparams)(((0.5.*concO1O1.*newparams(1).*(1 + ...
                newparams(6)./(8*concO1O1) - (1/8)*sqrt((newparams(6)./concO1O1).^2 + ...
                16*newparams(6)./concO1O1))./(newparams(4)^2))./(1+...
                2*concO1O1./newparams(4)+concO1O1.^2./(newparams(4)^2)+...
                0.5.*concO1O1.*newparams(1).*(1 + ...
                newparams(6)./(8*concO1O1) - (1/8)*sqrt((newparams(6)./concO1O1).^2 + ...
                16*newparams(6)./concO1O1))./(newparams(4)^2))-newmeansO1O1(j,:))./newSEsO1O1(j,:));
            newf_O2O1 = @(newparams)(((0.5.*concO2O1.*newparams(1).*(1 + ...
                newparams(6)./(8*concO2O1) - (1/8)*sqrt((newparams(6)./concO2O1).^2 + ...
                16*newparams(6)./concO2O1))./(newparams(4)*newparams(5)))./(1+...
                concO2O1./newparams(4)+concO2O1./newparams(5)+concO2O1.^2./(newparams(4)*newparams(5))+...
                0.5.*concO2O1.*newparams(1).*(1 + ...
                newparams(6)./(8*concO2O1) - (1/8)*sqrt((newparams(6)./concO2O1).^2 + ...
                16*newparams(6)./concO2O1))./(newparams(4)*newparams(5)))-newmeansO2O1(j,:))./newSEsO2O1(j,:));
            newf_TA = @(newparams)(((0.5.*concTA.*newparams(2).*(1 + ...
                newparams(6)./(8*concTA) - (1/8)*sqrt((newparams(6)./concTA).^2 + ...
                16*newparams(6)./concTA))./(newparams(3)*newparams(4)))./(1+...
                concTA./newparams(3)+concTA./newparams(4)+concTA.^2./(newparams(3)*newparams(4))+...
                0.5.*concTA.*newparams(2).*(1 + ...
                newparams(6)./(8*concTA) - (1/8)*sqrt((newparams(6)./concTA).^2 + ...
                16*newparams(6)./concTA))./(newparams(3)*newparams(4)))-newmeansTA(j,:))./newSEsTA(j,:));
        else
            newf_OidO1 = @(newparams)(((0.5.*concOidO1.*newparams(1)./(newparams(3)*newparams(4)))./(1+...
                concOidO1./newparams(3)+concOidO1./newparams(4)+concOidO1.^2./(newparams(3)*newparams(4))+...
                0.5.*concOidO1.*newparams(1)./(newparams(3)*newparams(4)))-newmeansOidO1(j,:))./newSEsOidO1(j,:));
            newf_O1O1 = @(newparams)(((0.5.*concO1O1.*newparams(1)./(newparams(4)^2))./(1+...
                2*concO1O1./newparams(4)+concO1O1.^2./(newparams(4)^2)+...
                0.5.*concO1O1.*newparams(1)./(newparams(4)^2))-newmeansO1O1(j,:))./newSEsO1O1(j,:));
            newf_O2O1 = @(newparams)(((0.5.*concO2O1.*newparams(1)./(newparams(4)*newparams(5)))./(1+...
                concO2O1./newparams(4)+concO2O1./newparams(5)+concO2O1.^2./(newparams(4)*newparams(5))+...
                0.5.*concO2O1.*newparams(1)./(newparams(4)*newparams(5)))-newmeansO2O1(j,:))./newSEsO2O1(j,:));
            newf_TA = @(newparams)(((0.5.*concTA.*newparams(2)./(newparams(3)*newparams(4)))./(1+...
                concTA./newparams(3)+concTA./newparams(4)+concTA.^2./(newparams(3)*newparams(4))+...
                0.5.*concTA.*newparams(2)./(newparams(3)*newparams(4)))-newmeansTA(j,:))./newSEsTA(j,:));
        end
        
        newmodel=@(newparams)([newf_OidO1(newparams) newf_O1O1(newparams) newf_O2O1(newparams) newf_TA(newparams)]);

        try
            if dimers
                paramsbs(j,:) = lsqnonlin(newmodel,[300,1000,6,40,200,0.1],[0 0 0 0 0 0],[],opt);
            else
                paramsbs(j,:) = lsqnonlin(newmodel,[300,1000,6,40,200],[0 0 0 0 0],[],opt);
            end
        catch
            missed = missed+1;
            oldparamsbs = paramsbs(1:j-1,:);
            clear paramsbs
            paramsbs=zeros(nboot-1,length(fitparams)); %I can't make the entries NaN's if the fitting fails because
                %std will be NaN then; so just remove that whole
                %attempt
            paramsbs(1:j-1,:) = oldparamsbs;
            keyboard
        end

        clear newmodel

%             if rem(j,10)==0
%                 disp(j)
%             end

    end
    paramSEs = std(paramsbs,1,1);
        
elseif strcmpi(dataset,'PUCBvsM')
    %params is [JB, JM, Kid, K1]
    paramnames{1} = 'JB';
    paramnames{2} = 'JM';
    paramnames{3} = 'Kid';
    paramnames{4} = 'K1';

    %same K:
    %params is [JB, JM, Kid, K1]
%         paramnames{1} = 'JB';
%         paramnames{2} = 'JM';
%         paramnames{3} = 'Kid';

    %enforce E8 Kd's:
%     paramnames{1} = 'JB';
%     paramnames{2} = 'JM';

    %Going to have the same problem here as in the for-loop, with the
    %bottom state in particular having a zero SE
    if length(find(SEsPUCB))~=length(SEsPUCB)
        SEsPUCB(~logical(SEsPUCB))=mean(SEsPUCB(logical(SEsPUCB)));
    end
    if length(find(SEsPUCM))~=length(SEsPUCM)
        SEsPUCM(~logical(SEsPUCM))=mean(SEsPUCM(logical(SEsPUCM)));
    end

    f_B = @(params)(((0.5.*concPUC.*params(1)./(params(3)*params(4)))./(1+...
        concPUC./params(3)+concPUC./params(4)+concPUC.^2./(params(3)*params(4))+...
        0.5.*concPUC.*params(1)./(params(3)*params(4))+...
        0.5.*concPUC.*params(2)./(params(3)*params(4)))-pLoopsPUCB)./SEsPUCB);
    f_M = @(params)(((0.5.*concPUC.*params(2)./(params(3)*params(4)))./(1+...
        concPUC./params(3)+concPUC./params(4)+concPUC.^2./(params(3)*params(4))+...
        0.5.*concPUC.*params(1)./(params(3)*params(4))+...
        0.5.*concPUC.*params(2)./(params(3)*params(4)))-pLoopsPUCM)./SEsPUCM);

    %same K:
%         f_B = @(params)(((0.5.*concPUC.*params(1)./(params(3)*params(3)))./(1+...
%             concPUC./params(3)+concPUC./params(3)+concPUC.^2./(params(3)*params(3))+...
%             0.5.*concPUC.*params(1)./(params(3)*params(3))+...
%             0.5.*concPUC.*params(2)./(params(3)*params(3)))-pLoopsPUCB)./SEsPUCB);
%         f_M = @(params)(((0.5.*concPUC.*params(2)./(params(3)*params(3)))./(1+...
%             concPUC./params(3)+concPUC./params(3)+concPUC.^2./(params(3)*params(3))+...
%             0.5.*concPUC.*params(1)./(params(3)*params(3))+...
%             0.5.*concPUC.*params(2)./(params(3)*params(3)))-pLoopsPUCM)./SEsPUCM);

    %enforce E8 Kd's:
%     f_B = @(params)(((0.5.*concPUC.*params(1)./(9.3*42.9))./(1+...
%         concPUC./9.3+concPUC./42.9+concPUC.^2./(9.3*42.9)+...
%         0.5.*concPUC.*params(1)./(9.3*42.9)+...
%         0.5.*concPUC.*params(2)./(9.3*42.9))-pLoopsPUCB)./SEsPUCB);
%     f_M = @(params)(((0.5.*concPUC.*params(2)./(9.3*42.9))./(1+...
%         concPUC./9.3+concPUC./42.9+concPUC.^2./(9.3*42.9)+...
%         0.5.*concPUC.*params(1)./(9.3*42.9)+...
%         0.5.*concPUC.*params(2)./(9.3*42.9))-pLoopsPUCM)./SEsPUCM);

    model=@(params)([f_B(params) f_M(params)]);

    %Use lsqnonlin to fit.  Second input are starting parameters, output is
    %best fit results for the params vector.  I can also add additional arguments
    %lower bound (vector), upper bound (vector) for all parameters.

    fitparams = lsqnonlin(model,[100,400,6,40],[0 0 0 0],[],opt);

    %same K:
    %fitparams = lsqnonlin(model,[100,400,6],[0 0 0],[],opt);

    %Enforce E8 Kd's:
    %fitparams = lsqnonlin(model,[100,400],[0 0],[],opt);


    %Find the errors:
    paramsbs = zeros(nboot,length(fitparams));
    for j=1:nboot

        %If any of the SEs are 0, that's a problem for the fitting
        %routine (div by 0!)  O2's SEs can be 0 for the lowest
        %concentration.  So, same as in fitconccurve: make the SE for
        %any zero points to be the average of the rest of the SEs:
        if length(find(newSEsPUCB(j,:)))~=length(newSEsPUCB(j,:))
            newSEsPUCB(j,~logical(newSEsPUCB(j,:)))=mean(newSEsPUCB(j,logical(newSEsPUCB(j,:))));
        end
        if length(find(newSEsPUCM(j,:)))~=length(newSEsPUCM(j,:))
            newSEsPUCM(j,~logical(newSEsPUCM(j,:)))=mean(newSEsPUCM(j,logical(newSEsPUCM(j,:))));
        end

        newf_PUCB = @(newparams)(((0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(4)))./(1+...
            concPUC./newparams(3)+concPUC./newparams(4)+concPUC.^2./(newparams(3)*newparams(4))+...
            0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(4))+...
            0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(4)))-newmeansPUCB(j,:))./newSEsPUCB(j,:));
        newf_PUCM = @(newparams)(((0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(4)))./(1+...
            concPUC./newparams(3)+concPUC./newparams(4)+concPUC.^2./(newparams(3)*newparams(4))+...
            0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(4))+...
            0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(4)))-newmeansPUCM(j,:))./newSEsPUCM(j,:));

        %Same K
%             newf_PUCB = @(newparams)(((0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(3)))./(1+...
%                 concPUC./newparams(3)+concPUC./newparams(3)+concPUC.^2./(newparams(3)*newparams(3))+...
%                 0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(3))+...
%                 0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(3)))-newmeansPUCB(j,:))./newSEsPUCB(j,:));
%             newf_PUCM = @(newparams)(((0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(3)))./(1+...
%                 concPUC./newparams(3)+concPUC./newparams(3)+concPUC.^2./(newparams(3)*newparams(3))+...
%                 0.5.*concPUC.*newparams(1)./(newparams(3)*newparams(3))+...
%                 0.5.*concPUC.*newparams(2)./(newparams(3)*newparams(3)))-newmeansPUCM(j,:))./newSEsPUCM(j,:));

%         newf_PUCB = @(newparams)(((0.5.*concPUC.*newparams(1)./(9.3*42.9))./(1+...
%             concPUC./9.3+concPUC./42.9+concPUC.^2./(9.3*42.9)+...
%             0.5.*concPUC.*newparams(1)./(9.3*42.9)+...
%             0.5.*concPUC.*newparams(2)./(9.3*42.9))-newmeansPUCB(j,:))./newSEsPUCB(j,:));
%         newf_PUCM = @(newparams)(((0.5.*concPUC.*newparams(2)./(9.3*42.9))./(1+...
%             concPUC./9.3+concPUC./42.9+concPUC.^2./(9.3*42.9)+...
%             0.5.*concPUC.*newparams(1)./(9.3*42.9)+...
%             0.5.*concPUC.*newparams(2)./(9.3*42.9))-newmeansPUCM(j,:))./newSEsPUCM(j,:));

        newmodel=@(newparams)([newf_PUCB(newparams) newf_PUCM(newparams)]);

        try
            paramsbs(j,:) = lsqnonlin(newmodel,[100,400,6,40],[0 0 0 0],[],opt);

            %same K
            %paramsbs(j,:) = lsqnonlin(newmodel,[100,400,6],[0 0 0],[],opt);

            %enforce E8 Kd's
            %paramsbs(j,:) = lsqnonlin(newmodel,[100,400],[0 0],[],opt);
        catch
            missed = missed+1;
            oldparamsbs = paramsbs(1:j-1,:);
            clear paramsbs
            paramsbs=zeros(nboot-1,length(fitparams)); %I can't make the entries NaN's if the fitting fails because
                %std will be NaN then; so just remove that whole
                %attempt
            paramsbs(1:j-1,:) = oldparamsbs;
            keyboard
        end

        clear newmodel

%             if rem(j,10)==0
%                 disp(j)
%             end

    end
    paramSEs = std(paramsbs,1,1);
    
elseif strcmpi(dataset,'E8OidO1O2TAPUC') || strcmpi(dataset,'E8OidO1O2TAPUCsub0') || strcmpi(dataset,'E8OidO1O2TAPUCsub02')
    %params is [JE8, JTA, JPUCB, JPUCM, Kid, K1, K2]
    paramnames{1} = 'JE8';
    paramnames{2} = 'JTA';
    paramnames{3} = 'JPUCB';
    paramnames{4} = 'JPUCM';
    paramnames{5} = 'Kid';
    paramnames{6} = 'K1';
    paramnames{7} = 'K2';

    %PUC data may have SE=0 for one of the loops
    if length(find(SEsPUCB))~=length(SEsPUCB)
        SEsPUCB(~logical(SEsPUCB))=mean(SEsPUCB(logical(SEsPUCB)));
    end
    if length(find(SEsPUCM))~=length(SEsPUCM)
        SEsPUCM(~logical(SEsPUCM))=mean(SEsPUCM(logical(SEsPUCM)));
    end
    
    if dimers
        paramnames{8} = 'KDT';
        
        f_OidO1 = @(params)(((0.5.*concOidO1.*params(1).*(1 + ...
            params(8)./(8*concOidO1) - (1/8)*sqrt((params(8)./concOidO1).^2 + ...
            16*params(8)./concOidO1))./(params(5)*params(6)))./(1+...
            concOidO1./params(5)+concOidO1./params(6)+concOidO1.^2./(params(5)*params(6))+...
            0.5.*concOidO1.*params(1).*(1 + ...
            params(8)./(8*concOidO1) - (1/8)*sqrt((params(8)./concOidO1).^2 + ...
            16*params(8)./concOidO1))./(params(5)*params(6)))-pLoopsOidO1)./SEsOidO1);
        f_O1O1 = @(params)(((0.5.*concO1O1.*params(1).*(1 + ...
            params(8)./(8*concO1O1) - (1/8)*sqrt((params(8)./concO1O1).^2 + ...
            16*params(8)./concO1O1))./(params(6)^2))./(1+...
            2*concO1O1./params(6)+concO1O1.^2./(params(6)^2)+...
            0.5.*concO1O1.*params(1).*(1 + ...
            params(8)./(8*concO1O1) - (1/8)*sqrt((params(8)./concO1O1).^2 + ...
            16*params(8)./concO1O1))./(params(6)^2))-pLoopsO1O1)./SEsO1O1);
       f_O2O1 = @(params)(((0.5.*concO2O1.*params(1).*(1 + ...
            params(8)./(8*concO2O1) - (1/8)*sqrt((params(8)./concO2O1).^2 + ...
            16*params(8)./concO2O1))./(params(6)*params(7)))./(1+...
            concO2O1./params(6)+concO2O1./params(7)+concO2O1.^2./(params(6)*params(7))+...
            0.5.*concO2O1.*params(1).*(1 + ...
            params(8)./(8*concO2O1) - (1/8)*sqrt((params(8)./concO2O1).^2 + ...
            16*params(8)./concO2O1))./(params(6)*params(7)))-pLoopsO2O1)./SEsO2O1); 
        f_TA = @(params)(((0.5.*concTA.*params(2).*(1 + ...
            params(8)./(8*concTA) - (1/8)*sqrt((params(8)./concTA).^2 + ...
            16*params(8)./concTA))./(params(5)*params(6)))./(1+...
            concTA./params(5)+concTA./params(6)+concTA.^2./(params(5)*params(6))+...
            0.5.*concTA.*params(2).*(1 + ...
            params(8)./(8*concTA) - (1/8)*sqrt((params(8)./concTA).^2 + ...
            16*params(8)./concTA))./(params(5)*params(6)))-pLoopsTA)./SEsTA);
        f_B = @(params)(((0.5.*concPUC.*params(3).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6)))./(1+...
            concPUC./params(5)+concPUC./params(6)+concPUC.^2./(params(5)*params(6))+...
            0.5.*concPUC.*params(3).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6))+...
            0.5.*concPUC.*params(4).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6)))-pLoopsPUCB)./SEsPUCB);
        f_M = @(params)(((0.5.*concPUC.*params(4).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6)))./(1+...
            concPUC./params(5)+concPUC./params(6)+concPUC.^2./(params(5)*params(6))+...
            0.5.*concPUC.*params(3).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6))+...
            0.5.*concPUC.*params(4).*(1 + ...
            params(8)./(8*concPUC) - (1/8)*sqrt((params(8)./concPUC).^2 + ...
            16*params(8)./concPUC))./(params(5)*params(6)))-pLoopsPUCM)./SEsPUCM);
    else
        f_OidO1 = @(params)(((0.5.*concOidO1.*params(1)./(params(5)*params(6)))./(1+...
            concOidO1./params(5)+concOidO1./params(6)+concOidO1.^2./(params(5)*params(6))+...
            0.5.*concOidO1.*params(1)./(params(5)*params(6)))-pLoopsOidO1)./SEsOidO1);
        f_O1O1 = @(params)(((0.5.*concO1O1.*params(1)./(params(6)^2))./(1+...
            2*concO1O1./params(6)+concO1O1.^2./(params(6)^2)+...
            0.5.*concO1O1.*params(1)./(params(6)^2))-pLoopsO1O1)./SEsO1O1);
        f_O2O1 = @(params)(((0.5.*concO2O1.*params(1)./(params(6)*params(7)))./(1+...
            concO2O1./params(6)+concO2O1./params(7)+concO2O1.^2./(params(6)*params(7))+...
            0.5.*concO2O1.*params(1)./(params(6)*params(7)))-pLoopsO2O1)./SEsO2O1);
        f_TA = @(params)(((0.5.*concTA.*params(2)./(params(5)*params(6)))./(1+...
            concTA./params(5)+concTA./params(6)+concTA.^2./(params(5)*params(6))+...
            0.5.*concTA.*params(2)./(params(5)*params(6)))-pLoopsTA)./SEsTA);
        f_B = @(params)(((0.5.*concPUC.*params(3)./(params(5)*params(6)))./(1+...
            concPUC./params(5)+concPUC./params(6)+concPUC.^2./(params(5)*params(6))+...
            0.5.*concPUC.*params(3)./(params(5)*params(6))+...
            0.5.*concPUC.*params(4)./(params(5)*params(6)))-pLoopsPUCB)./SEsPUCB);
        f_M = @(params)(((0.5.*concPUC.*params(4)./(params(5)*params(6)))./(1+...
            concPUC./params(5)+concPUC./params(6)+concPUC.^2./(params(5)*params(6))+...
            0.5.*concPUC.*params(3)./(params(5)*params(6))+...
            0.5.*concPUC.*params(4)./(params(5)*params(6)))-pLoopsPUCM)./SEsPUCM);
    end
    
    model=@(params)([f_OidO1(params) f_O1O1(params) f_O2O1(params) f_TA(params) f_B(params) f_M(params)]);

    if dimers
        fitparams = lsqnonlin(model,[300,1000,200,400,6,40,200,0.1],[0 0 0 0 0 0 0 0],[],opt);
    else
        fitparams = lsqnonlin(model,[300,1000,200,400,6,40,200],[0 0 0 0 0 0 0],[],opt);
    end

    %Find the errors:
    paramsbs = zeros(nboot,length(fitparams));
    for j=1:nboot

        %If any of the SEs are 0, that's a problem for the fitting
        %routine (div by 0!)  O2's SEs can be 0 for the lowest
        %concentration.  So, same as in fitconccurve: make the SE for
        %any zero points to be the average of the rest of the SEs:
        if length(find(newSEsOidO1(j,:)))~=length(newSEsOidO1(j,:))
            newSEsOidO1(j,~logical(newSEsOidO1(j,:)))=mean(newSEsOidO1(j,logical(newSEsOidO1(j,:))));
        end
        if length(find(newSEsO1O1(j,:)))~=length(newSEsO1O1(j,:))
            newSEsO1O1(j,~logical(newSEsO1O1(j,:)))=mean(newSEsO1O1(j,logical(newSEsO1O1(j,:))));
        end
        if length(find(newSEsO2O1(j,:)))~=length(newSEsO2O1(j,:))
            newSEsO2O1(j,~logical(newSEsO2O1(j,:)))=mean(newSEsO2O1(j,logical(newSEsO2O1(j,:))));
        end
        if length(find(newSEsTA(j,:)))~=length(newSEsTA(j,:))
            newSEsTA(j,~logical(newSEsTA(j,:)))=mean(newSEsTA(j,logical(newSEsTA(j,:))));
        end
        if length(find(newSEsPUCB(j,:)))~=length(newSEsPUCB(j,:))
            newSEsPUCB(j,~logical(newSEsPUCB(j,:)))=mean(newSEsPUCB(j,logical(newSEsPUCB(j,:))));
        end
        if length(find(newSEsPUCM(j,:)))~=length(newSEsPUCM(j,:))
            newSEsPUCM(j,~logical(newSEsPUCM(j,:)))=mean(newSEsPUCM(j,logical(newSEsPUCM(j,:))));
        end

        if dimers
            newf_OidO1 = @(newparams)(((0.5.*concOidO1.*newparams(1).*(1 + ...
                newparams(8)./(8*concOidO1) - (1/8)*sqrt((newparams(8)./concOidO1).^2 + ...
                16*newparams(8)./concOidO1))./(newparams(5)*newparams(6)))./(1+...
                concOidO1./newparams(5)+concOidO1./newparams(6)+concOidO1.^2./(newparams(5)*newparams(6))+...
                0.5.*concOidO1.*newparams(1).*(1 + ...
                newparams(8)./(8*concOidO1) - (1/8)*sqrt((newparams(8)./concOidO1).^2 + ...
                16*newparams(8)./concOidO1))./(newparams(5)*newparams(6)))-newmeansOidO1(j,:))./newSEsOidO1(j,:));
            newf_O1O1 = @(newparams)(((0.5.*concO1O1.*newparams(1).*(1 + ...
                newparams(8)./(8*concO1O1) - (1/8)*sqrt((newparams(8)./concO1O1).^2 + ...
                16*newparams(8)./concO1O1))./(newparams(6)^2))./(1+...
                2*concO1O1./newparams(6)+concO1O1.^2./(newparams(6)^2)+...
                0.5.*concO1O1.*newparams(1).*(1 + ...
                newparams(8)./(8*concO1O1) - (1/8)*sqrt((newparams(8)./concO1O1).^2 + ...
                16*newparams(8)./concO1O1))./(newparams(6)^2))-newmeansO1O1(j,:))./newSEsO1O1(j,:));
            newf_O2O1 = @(newparams)(((0.5.*concO2O1.*newparams(1).*(1 + ...
                newparams(8)./(8*concO2O1) - (1/8)*sqrt((newparams(8)./concO2O1).^2 + ...
                16*newparams(8)./concO2O1))./(newparams(6)*newparams(7)))./(1+...
                concO2O1./newparams(6)+concO2O1./newparams(7)+concO2O1.^2./(newparams(6)*newparams(7))+...
                0.5.*concO2O1.*newparams(1).*(1 + ...
                newparams(8)./(8*concO2O1) - (1/8)*sqrt((newparams(8)./concO2O1).^2 + ...
                16*newparams(8)./concO2O1))./(newparams(6)*newparams(7)))-newmeansO2O1(j,:))./newSEsO2O1(j,:)); 
            newf_TA = @(newparams)(((0.5.*concTA.*newparams(2).*(1 + ...
                newparams(8)./(8*concTA) - (1/8)*sqrt((newparams(8)./concTA).^2 + ...
                16*newparams(8)./concTA))./(newparams(5)*newparams(6)))./(1+...
                concTA./newparams(5)+concTA./newparams(6)+concTA.^2./(newparams(5)*newparams(6))+...
                0.5.*concTA.*newparams(2).*(1 + ...
                newparams(8)./(8*concTA) - (1/8)*sqrt((newparams(8)./concTA).^2 + ...
                16*newparams(8)./concTA))./(newparams(5)*newparams(6)))-newmeansTA(j,:))./newSEsTA(j,:));
            newf_B = @(newparams)(((0.5.*concPUC.*newparams(3).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6)))./(1+...
                concPUC./newparams(5)+concPUC./newparams(6)+concPUC.^2./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(3).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(4).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6)))-newmeansPUCB(j,:))./newSEsPUCB(j,:));
            newf_M = @(newparams)(((0.5.*concPUC.*newparams(4).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6)))./(1+...
                concPUC./newparams(5)+concPUC./newparams(6)+concPUC.^2./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(3).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(4).*(1 + ...
                newparams(8)./(8*concPUC) - (1/8)*sqrt((newparams(8)./concPUC).^2 + ...
                16*newparams(8)./concPUC))./(newparams(5)*newparams(6)))-newmeansPUCM(j,:))./newSEsPUCM(j,:));
        else
            newf_OidO1 = @(newparams)(((0.5.*concOidO1.*newparams(1)./(newparams(5)*newparams(6)))./(1+...
                concOidO1./newparams(5)+concOidO1./newparams(6)+concOidO1.^2./(newparams(5)*newparams(6))+...
                0.5.*concOidO1.*newparams(1)./(newparams(5)*newparams(6)))-newmeansOidO1(j,:))./newSEsOidO1(j,:));
            newf_O1O1 = @(newparams)(((0.5.*concO1O1.*newparams(1)./(newparams(6)^2))./(1+...
                2*concO1O1./newparams(6)+concO1O1.^2./(newparams(6)^2)+...
                0.5.*concO1O1.*newparams(1)./(newparams(6)^2))-newmeansO1O1(j,:))./newSEsO1O1(j,:));
            newf_O2O1 = @(newparams)(((0.5.*concO2O1.*newparams(1)./(newparams(6)*newparams(7)))./(1+...
                concO2O1./newparams(6)+concO2O1./newparams(7)+concO2O1.^2./(newparams(6)*newparams(7))+...
                0.5.*concO2O1.*newparams(1)./(newparams(6)*newparams(7)))-newmeansO2O1(j,:))./newSEsO2O1(j,:));
            newf_TA = @(newparams)(((0.5.*concTA.*newparams(2)./(newparams(5)*newparams(6)))./(1+...
                concTA./newparams(5)+concTA./newparams(6)+concTA.^2./(newparams(5)*newparams(6))+...
                0.5.*concTA.*newparams(2)./(newparams(5)*newparams(6)))-newmeansTA(j,:))./newSEsTA(j,:));
            newf_B = @(newparams)(((0.5.*concPUC.*newparams(3)./(newparams(5)*newparams(6)))./(1+...
                concPUC./newparams(5)+concPUC./newparams(6)+concPUC.^2./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(3)./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(4)./(newparams(5)*newparams(6)))-newmeansPUCB(j,:))./newSEsPUCB(j,:));
            newf_M = @(newparams)(((0.5.*concPUC.*newparams(4)./(newparams(5)*newparams(6)))./(1+...
                concPUC./newparams(5)+concPUC./newparams(6)+concPUC.^2./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(3)./(newparams(5)*newparams(6))+...
                0.5.*concPUC.*newparams(4)./(newparams(5)*newparams(6)))-newmeansPUCM(j,:))./newSEsPUCM(j,:));
        end
        
        newmodel=@(newparams)([newf_OidO1(newparams) newf_O1O1(newparams) newf_O2O1(newparams) newf_TA(newparams) newf_B(newparams) newf_M(newparams)]);

        try
            if dimers
                paramsbs(j,:) = lsqnonlin(newmodel,[300,1000,200,400,6,40,200,0.1],[0 0 0 0 0 0 0 0],[],opt);
            else
                paramsbs(j,:) = lsqnonlin(newmodel,[300,1000,200,400,6,40,200],[0 0 0 0 0 0 0],[],opt);
            end
        catch
            missed = missed+1;
            oldparamsbs = paramsbs(1:j-1,:);
            clear paramsbs
            paramsbs=zeros(nboot-1,length(fitparams)); %I can't make the entries NaN's if the fitting fails because
                %std will be NaN then; so just remove that whole
                %attempt
            paramsbs(1:j-1,:) = oldparamsbs;
            keyboard
        end

        clear newmodel

%             if rem(j,10)==0
%                 disp(j)
%             end

    end
    paramSEs = std(paramsbs,1,1);

elseif strcmpi(dataset,'E8107BvsMsub02')

    paramnames{1} = 'JB';
    paramnames{2} = 'JM';

    %Going to have the same problem here as in the for-loop, with the
    %bottom state in particular having a zero SE
    if length(find(SEsE107B))~=length(SEsE107B)
        SEsE107B(~logical(SEsE107B))=mean(SEsE107B(logical(SEsE107B)));
    end
    if length(find(SEsE107M))~=length(SEsE107M)
        SEsE107M(~logical(SEsE107M))=mean(SEsE107M(logical(SEsE107M)));
    end

    %Being lazy since concE107 = concE107B = concE107M
    f_B = @(params)(((0.5.*concE107.*params(1)./(Kid*K1))./(1+...
        concE107./Kid+concE107./K1+concE107.^2./(Kid*K1)+...
        0.5.*concE107.*params(1)./(Kid*K1)+...
        0.5.*concE107.*params(2)./(Kid*K1))-pLoopsE107B)./SEsE107B);
    f_M = @(params)(((0.5.*concE107.*params(2)./(Kid*K1))./(1+...
        concE107./Kid+concE107./K1+concE107.^2./(Kid*K1)+...
        0.5.*concE107.*params(1)./(Kid*K1)+...
        0.5.*concE107.*params(2)./(Kid*K1))-pLoopsE107M)./SEsE107M);

    model=@(params)([f_B(params) f_M(params)]);

    %Use lsqnonlin to fit.  Second input are starting parameters, output is
    %best fit results for the params vector.  I can also add additional arguments
    %lower bound (vector), upper bound (vector) for all parameters.

    fitparams = lsqnonlin(model,[100,200],[0 0],[],opt);


    %Find the errors:
    paramsbs = zeros(nboot,length(fitparams));
    for j=1:nboot

        %If any of the SEs are 0, that's a problem for the fitting
        %routine (div by 0!)  O2's SEs can be 0 for the lowest
        %concentration.  So, same as in fitconccurve: make the SE for
        %any zero points to be the average of the rest of the SEs:
        if length(find(newSEsE107B(j,:)))~=length(newSEsE107B(j,:))
            newSEsE107B(j,~logical(newSEsE107B(j,:)))=mean(newSEsE107B(j,logical(newSEsE107B(j,:))));
        end
        if length(find(newSEsE107M(j,:)))~=length(newSEsE107M(j,:))
            newSEsE107M(j,~logical(newSEsE107M(j,:)))=mean(newSEsE107M(j,logical(newSEsE107M(j,:))));
        end

        newf_E107B = @(newparams)(((0.5.*concE107.*newparams(1)./(Kiddistrib(j)*K1distrib(j)))./(1+...
            concE107./Kiddistrib(j)+concE107./K1distrib(j)+concE107.^2./(Kiddistrib(j)*K1distrib(j))+...
            0.5.*concE107.*newparams(1)./(Kiddistrib(j)*K1distrib(j))+...
            0.5.*concE107.*newparams(2)./(Kiddistrib(j)*Kiddistrib(j)))-newmeansE107B(j,:))./newSEsE107B(j,:));
        newf_E107M = @(newparams)(((0.5.*concE107.*newparams(2)./(Kiddistrib(j)*K1distrib(j)))./(1+...
            concE107./Kiddistrib(j)+concE107./K1distrib(j)+concE107.^2./(Kiddistrib(j)*K1distrib(j))+...
            0.5.*concE107.*newparams(1)./(Kiddistrib(j)*K1distrib(j))+...
            0.5.*concE107.*newparams(2)./(Kiddistrib(j)*K1distrib(j)))-newmeansE107M(j,:))./newSEsE107M(j,:));

        newmodel=@(newparams)([newf_E107B(newparams) newf_E107M(newparams)]);

        try
            paramsbs(j,:) = lsqnonlin(newmodel,[100,200],[0 0],[],opt);
        catch
            missed = missed+1;
            oldparamsbs = paramsbs(1:j-1,:);
            clear paramsbs
            paramsbs=zeros(nboot-1,length(fitparams)); %I can't make the entries NaN's if the fitting fails because
                %std will be NaN then; so just remove that whole
                %attempt
            paramsbs(1:j-1,:) = oldparamsbs;
            keyboard
        end

        clear newmodel

%             if rem(j,10)==0
%                 disp(j)
%             end

    end
    paramSEs = std(paramsbs,1,1);

    
end

disp(missed)

%Figures to check bootstrapping+fitting part:
xvals = 50:50:floor(nboot/50)*50;
yvals = zeros(length(fitparams),length(xvals));

for k=1:length(xvals)
    yvals(:,k) = std(paramsbs(1:k*50,:),1,1);
end

for p=1:length(fitparams)              
    figure
    [n,xout] = hist(paramsbs(:,p),nbins);
    bar(xout,n)
    title(strcat(paramnames{p},'Distribution'))

    figure
    plot(xvals,yvals(p,:),'-x')
    title(strcat('Error on ',paramnames{p}))
    xlabel('Number of Resamplings')

    clear n xout
end

if exist('/Users/Steph/Desktop/bootstrapping')
    save('/Users/Steph/Desktop/bootstrapping2.mat','paramsbs','paramnames','fitparams','paramSEs')
else
    save('/Users/Steph/Desktop/bootstrapping.mat','paramsbs','paramnames','fitparams','paramSEs')
end
    
    
%Subroutine that bootstraps the data
    function [newmeans,newSEs,totmeans,totSEs] = bootstrap_concs(conc,data)
    
        newmeans = zeros(nboot,length(conc)); %New means will have a column for
            %each concentration, where each column is the output of
            %bootstrp in the next for loop
        newSEs = zeros(nboot,length(conc));
        totmeans = zeros(1,length(conc));
        totSEs = zeros(1,length(conc),1);
        newstats=cell(length(conc),1);
        newdistribs=cell(length(conc),1);

        for i=1:length(conc) %For each concentation
            totmeans(i) = mean(data{i}); %First get the mean of the total distribution
            totSEs(i) = std(data{i})/sqrt(length(data{i})-1); %And the SE, for weights for the fits
            [newstats{i},newdistribs{i}] = bootstrp(nboot,@(x)[mean(x), std(x)/sqrt(length(x)-1)],data{i}'); %data{i} needs to be a column vector
            newmeans(:,i) = newstats{i}(:,1);
            newSEs(:,i) = newstats{i}(:,2);
        end
    end
end

