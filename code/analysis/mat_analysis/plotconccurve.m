%function [h,PlotHandle] = plotconccurve(concs, pLoops, SEs, legendlabels, theoryx,theoryy, varargin)
%
%For plotting concentration titrations, since this comes up so much.
%Returns a handle to the figure and the plot elements.
%
%Inputs:
% concs: a cell array of data, each element one x-set
% pLoops: a cell array same length as concs
% SEs: cell array same length as concs
% legend: user-inputed legend entries; must be a cell array of strings
% theoryx, theoryy: these are R-vs-pLoops to plot as a curve (not data points with
%   errors); again a cell array.  need not be the same length as concs.
%   If no theory to plot, enter [] for both.
% Optional input: color and line specs for data and/or theory.  If want to
% specify theory but not data, enter [] for the first optional input.
% Options are:
%	'E8",'TA','PUC','PUCM','PUCB': these use default colors and line specs for these DNA types
%   {'xb','.r',...}: cell array same length as concs or theory with
%       user-specified input for each x,y pair.  THESE MUST HAVE THE COLOR
%       LAST!
% second optional input added 3/2011: 'Lengths' indicates this is a length
%   series, not a concentration curve, and adjusts the xlabel etc
%   accordingly.  Similarly with 'lengthsHG','JHG' and 'FHG'.  
%   If no theory color specs, enter [], then '<input>'.  You
%   must also specify a third optional input:
% third optional input added 9/2011: Add 'OpDist' to use operator 
%   center-to-center spacing; else 'LoopLen'for just the length of the loop, 
%   no operators.
%
%Steph 10/10

function [h,PlotHandle] = plotconccurve(concs, pLoops, SEs, legendlabels, theoryx,theoryy, varargin)

%Determine the linespecs for the data and the theory:

if nargin>6
    if ~isempty(varargin{1}) %Means user had opinions about datastyle
        datastyle = varargin{1};
        for i=1:length(datastyle)
            if strcmpi(datastyle{i},'E8')
                datastyle{i} = '.k';
            elseif strcmpi(datastyle{i},'TA')
                datastyle{i} = '.r';
            elseif strcmpi(datastyle{i},'PUC')
                datastyle{i} = '.g';
            elseif strcmpi(datastyle{i},'PUCM');
                datastyle{i} = 'xg';
            elseif strcmpi(datastyle{i},'PUCB')
                datastyle{i} = 'og';
            elseif strcmpi(datastyle{i},'E8O1O1')
                datastyle{i} = 'xb';
            elseif strcmpi(datastyle{i},'E8O2O1')
                datastyle{i} = 'om';
            end
        end
    end
        
    if length(varargin)>1 %Means user  had opinions about theorystyle also
        theorystyle = varargin{2};
        for i=1:length(theorystyle)
            if strcmpi(theorystyle{i},'E8')
                theorystyle{i} = '-k';
            elseif strcmpi(theorystyle{i},'TA')
                theorystyle{i} = '-r';
            elseif strcmpi(theorystyle{i},'PUC')
                theorystyle{i} = '-g';
            elseif strcmpi(theorystyle{i},'PUCM');
                theorystyle{i} = ':g';
            elseif strcmpi(theorystyle{i},'PUCB')
                theorystyle{i} = '--g';
            elseif strcmpi(theorystyle{i},'E8O1O1')
                theorystyle{i} = '-b';
            elseif strcmpi(theorystyle{i},'E8O2O1')
                theorystyle{i} = '-m';
            end
        end
    else %Means user only cared about datastyle, or there isn't any theory to plot
        if ~isempty(theoryx) %Not only data points
            theorystyle{1} = '-k';
            theorystyle{2} = '-r';
            theorystyle{3} = '-g';
            theorystyle{4} = '-b';
            theorystyle{5} = '-y';
        end
    end
else %User doesn't care about data or theory style
    datastyle{1} = '.k';
    datastyle{2} = '.r';
    datastyle{3} = '.g';
    datastyle{4} = '.b';
    datastyle{5} = '.y'; %If I plot more than 5 things at once no one will be able to decipher it!
    if ~isempty(theoryx) %Not only data points
        theorystyle{1} = '-k';
        theorystyle{2} = '-r';
        theorystyle{3} = '-g';
        theorystyle{4} = '-b';
        theorystyle{5} = '-y';
    end
end
        

h=figure;

%Plot just concs vs ploops to get the legend right
PlotHandle = plot(concs{1},pLoops{1},datastyle{1});
hold on

for i=2:length(concs)
    PlotHandle(end+1) = plot(concs{i},pLoops{i},datastyle{i});
end

%Plot any theory curves
if ~isempty(theoryx)
    PlotHandle(end+1) = plot(theoryx{1},theoryy{1},theorystyle{1});
    for i=2:length(theoryx)
        PlotHandle(end+1) = plot(theoryx{i},theoryy{i},theorystyle{i});
    end
end

%Replot the data, this time with errorbars
PlotHandle(end+1:end+length(concs{1})+1)=errorbarxyHG(concs{1},pLoops{1},[],...
    SEs{1},[],[],datastyle{1},datastyle{1}(end));

for i=2:length(concs)
    PlotHandle(end+1:end+length(concs{i})+1)=errorbarxyHG(concs{i},pLoops{i},[],...
        SEs{i},[],[],datastyle{i},datastyle{i}(end));
end

%Take care of the rest of the plot
hold off

if nargin>8 && strcmpi(varargin{3},'lengths')
    
    if strcmpi(varargin{4},'OpDist')
        xlabel('Distance between operators (bp)','FontSize',16)
        xlim([88+20.5 117+20.5])
    else
        xlabel('Loop Length (bp)','FontSize',16)
        xlim([88 117])
    end
    
    ylabel('Looping Probability','FontSize',16)
    ylim([0 1])
    
elseif nargin>8 && strcmpi(varargin{3},'lengthsHG')
    
    if strcmpi(varargin{4},'OpDist')
        xlabel('Distance between operators (bp)','FontSize',16)
        xlim([55+36+20.5 89+36+20.5])
    else
        xlabel('Loop Length (bp)','FontSize',16)
        xlim([55+36 89+36])
    end

    ylabel('Looping Probability','FontSize',16)
    ylim([0 1])
    
elseif nargin>8 && strcmpi(varargin{3},'JHG')
    if strcmpi(varargin{4},'OpDist')
        xlabel('Distance between operators (bp)','FontSize',16)
        xlim([55+36+20.5 89+36+20.5])
    else
        xlabel('Loop Length (bp)','FontSize',16)
        xlim([55+36 89+36])
    end
    
    ylabel('J_{loop} (M)','FontSize',16)
    set(gca,'YScale','log')
    
elseif nargin>8 && strcmpi(varargin{3},'FHG')
    
    if strcmpi(varargin{4},'OpDist')
        xlabel('Distance between operators (bp)','FontSize',16)
        xlim([55+36+20.5 89+36+20.5])
    else
        xlabel('Loop Length (bp)','FontSize',16)
        xlim([55+36 89+36])
    end
    
    ylabel('\Delta F_{loop} (kT)','FontSize',16)

else
    set(gca,'XScale','log')
    xlabel('LacI Concentration (M)','FontSize',16)
    xlim([10^-15 10^-5])
    set(gca,'XTick',[10^-15 10^-14 10^-13 10^-12 10^-11 10^-10 10^-9 10^-8 10^-7 10^-6 10^-5])
    ylabel('Looping Probability','FontSize',16)
    ylim([0 1])
end

set(gca,'FontSize',14)
legend(legendlabels)

StandardFigure(PlotHandle,gca)