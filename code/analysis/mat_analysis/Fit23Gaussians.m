%Steph 11/08

%Modified version of HGFit23Gaussians to interface with
%singlebeadanalysisV2.

function GaussFit=Fit23Gaussians(newxyr,Name,startp,Identity,AllFit,xx)

%Identity parameters:
%B: Bottom loop
%M: Middle loop
%U: Unlooped state
Identities={'B','M','U'};

if isempty(AllFit)
    AllFit.b1=startp(2);
    AllFit.b2=startp(5);
    AllFit.b3=startp(8);
end
    

%xx=[100:2:260];
%xxRange=linspace(100,260);

%For Wtlac:
% xx=[80:2:240];
% xxRange=linspace(80,240);

%GaussFit.newr=newxyr;

[n,xx]=hist(newxyr,xx); %xx are the bin locations, n is the distribution
all=n./sum(n(2:end));

GaussFit.n=n;
GaussFit.all=all;

GaussFit.Approved=1;
GaussFit.name=Name;
GaussFit.startp=startp;



if length(startp)==9
    %opts = fitoptions('gauss3','start',startp,'Algorithm','Trust-Region');
    %Steph 1/18/12: Changed the above lien to the following so that it
    %won't fit negative parameter values
    opts = fitoptions('gauss3','start',startp,'Algorithm','Trust-Region','lower',[0 0 0]);
    GaussFit.fp = fit(xx(2:end)',all(2:end)','gauss3',opts);
    
    %Calculates the looping probability by finding the minimum between
    %Gaussians and integrating the fitted functions.
    GaussParameters=[GaussFit.fp.a1,GaussFit.fp.a2,GaussFit.fp.a3;...
        GaussFit.fp.b1,GaussFit.fp.b2,GaussFit.fp.b3;...
        GaussFit.fp.c1,GaussFit.fp.c2,GaussFit.fp.c3];

    %Steph 9/10: THIS IS A BUG! it potentially scrambles p relative to the
    %peak identities and screws up the pLoop calculated for this bead!
    %[GaussPositions,IX]=sort(GaussParameters(2,:));

    %Integrate the area under each Gaussian
    syms x
    %Steph changed this 9/10: this scrambles the connection between p,
    %Identity, and the fit parameters for each peak
%     GaussSymbolic=GaussParameters(1,IX(1))*...
%         exp(-((x-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2)+...
%         GaussParameters(1,IX(2))*...
%         exp(-((x-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2)+...
%         GaussParameters(1,IX(3))*...
%         exp(-((x-GaussParameters(2,IX(3)))/GaussParameters(3,IX(3))).^2);
%     Gauss1Symbolic=GaussParameters(1,IX(1))*...
%         exp(-((x-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2);
%     Gauss2Symbolic=GaussParameters(1,IX(2))*...
%         exp(-((x-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2);
%     Gauss3Symbolic=GaussParameters(1,IX(3))*...
%         exp(-((x-GaussParameters(2,IX(3)))/GaussParameters(3,IX(3))).^2);

    GaussSymbolic=GaussParameters(1,1)*...
        exp(-((x-GaussParameters(2,1))/GaussParameters(3,1)).^2)+...
        GaussParameters(1,2)*...
        exp(-((x-GaussParameters(2,2))/GaussParameters(3,2)).^2)+...
        GaussParameters(1,3)*...
        exp(-((x-GaussParameters(2,3))/GaussParameters(3,3)).^2);
    Gauss1Symbolic=GaussParameters(1,1)*...
        exp(-((x-GaussParameters(2,1))/GaussParameters(3,1)).^2);
    Gauss2Symbolic=GaussParameters(1,2)*...
        exp(-((x-GaussParameters(2,2))/GaussParameters(3,2)).^2);
    Gauss3Symbolic=GaussParameters(1,3)*...
        exp(-((x-GaussParameters(2,3))/GaussParameters(3,3)).^2);

    GausstAll=eval(int(GaussSymbolic,x,xx(1),xx(end)));
    GaussIndept1=eval(int(Gauss1Symbolic,x,xx(1),xx(end)));
    GaussIndept2=eval(int(Gauss2Symbolic,x,xx(1),xx(end)));
    GaussIndept3=eval(int(Gauss3Symbolic,x,xx(1),xx(end)));
    GaussIndept23=GaussIndept2+GaussIndept3;

    %GaussFit.ploop=GaussIndept23/GausstAll;
    p1=GaussIndept1/GausstAll;
    p2=GaussIndept2/GausstAll;
    p3=GaussIndept3/GausstAll;
    
    GaussFit.p=[p1,p2,p3];
    
    if isempty(Identity)
        for i=1:length(GaussFit.p)
            eval(['[Min,Pos]=min([abs(AllFit.b1-GaussFit.fp.b',num2str(i),...
                '),abs(AllFit.b2-GaussFit.fp.b',num2str(i),...
                '),abs(AllFit.b3-GaussFit.fp.b',num2str(i),')]);']);
            Identity{i}=Identities{Pos};
        end
    end
    
elseif length(startp)==6
    %opts = fitoptions('gauss2','start',startp,'Algorithm','Trust-Region');
    %Steph 1/18/12: Changed the above lien to the following so that it
    %won't fit negative parameter values
    opts = fitoptions('gauss2','start',startp,'Algorithm','Trust-Region','lower',[0 0]);
    GaussFit.fp = fit(xx(2:end)',all(2:end)','gauss2',opts);
    
    %Calculates the looping probability by finding the minimum between
    %Gaussians and integrating the fitted functions.
    GaussParameters=[GaussFit.fp.a1,GaussFit.fp.a2;...
        GaussFit.fp.b1,GaussFit.fp.b2;...
        GaussFit.fp.c1,GaussFit.fp.c2];

    %Steph 9/10: THIS IS A BUG! it potentially scrambles p relative to the
    %peak identities and screws up the pLoop calculated for this bead!
    %[GaussPositions,IX]=sort(GaussParameters(2,:));

    %Integrate the area under each Gaussian
    syms x
%     GaussSymbolic=GaussParameters(1,IX(1))*...
%         exp(-((x-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2)+...
%         GaussParameters(1,IX(2))*...
%         exp(-((x-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2);
%     Gauss1Symbolic=GaussParameters(1,IX(1))*...
%         exp(-((x-GaussParameters(2,IX(1)))/GaussParameters(3,IX(1))).^2);
%     Gauss2Symbolic=GaussParameters(1,IX(2))*...
%         exp(-((x-GaussParameters(2,IX(2)))/GaussParameters(3,IX(2))).^2);
    
    GaussSymbolic=GaussParameters(1,1)*...
        exp(-((x-GaussParameters(2,1))/GaussParameters(3,1)).^2)+...
        GaussParameters(1,2)*...
        exp(-((x-GaussParameters(2,2))/GaussParameters(3,2)).^2);
    Gauss1Symbolic=GaussParameters(1,1)*...
        exp(-((x-GaussParameters(2,1))/GaussParameters(3,1)).^2);
    Gauss2Symbolic=GaussParameters(1,2)*...
        exp(-((x-GaussParameters(2,2))/GaussParameters(3,2)).^2);

    GausstAll=eval(int(GaussSymbolic,x,xx(1),xx(end)));
    GaussIndept1=eval(int(Gauss1Symbolic,x,xx(1),xx(end)));
    GaussIndept2=eval(int(Gauss2Symbolic,x,xx(1),xx(end)));


    %GaussFit.ploop=GaussIndept23/GausstAll;
    p1=GaussIndept1/GausstAll;
    p2=GaussIndept2/GausstAll;
    
    GaussFit.p=[p1,p2];
    
    if isempty(Identity)
        for i=1:length(GaussFit.p)
            eval(['[Min,Pos]=min([abs(AllFit.b1-GaussFit.fp.b',num2str(i),...
                '),abs(AllFit.b2-GaussFit.fp.b',num2str(i),')]);']);
            Identity{i}=Identities{Pos};
        end
    end
end

GaussFit.Identity=Identity;