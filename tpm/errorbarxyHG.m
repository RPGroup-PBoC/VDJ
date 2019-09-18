function PlotHandle=errorbarxyHG(x,y,lx,ly,ux,uy,linecol,errorcol,varargin)
%This function allows the user to plot the graph of x against y, along with both x and y errorbars.
%For the x and y errors it is possible to input both lower (lx and ly)  and upper  (ux and uy) values for the
%errors at a particular point.  If the upper values are not specified then the program assumes the errors 
%are symmetrical and use the lower values.  it is also possible to specify the plot line colour, marker, and 
%linestyle using the standard 'plot' command notation in the input variable 'linecol'.  Also the line colour for 
%the errobars can be specified in the variable 'errorcol'.  It is important to note that if these colour options 
%are to be used and any of the error limit vectors are empty then they should not be excluded, but presented 
%in a [] form signifying an empty vector.
%
%James Rooney,  17 October 2003

%2009/06/21: Modified by Hernan Garcia to output a handle.
%2009/06/28: Modified by Hernan Garcia to not show the perpendicular lines
%at the end of the errorbar. Look for (*) where commented out.
%2009/07/09: Modified such that more properties of the line can be input.

if exist('linecol','var')==0 | isempty(linecol)
    linecol='b';
end


if exist('errorcol','var')==0 | isempty(errorcol)
    errorcol='r';
end

if length(varargin)==2
    PlotHandle=plot(x,y,linecol,varargin{1},varargin{2});
else
    PlotHandle=plot(x,y,linecol);
end

if ishold
    KeepHold=1;
else
    KeepHold=0;
    hold on
end

xw=(max(x)-min(x))/100;
yw=(max(y)-min(y))/100;


lye=exist('ly','var');
lxe=exist('lx','var');
uye=exist('uy','var');
uxe=exist('ux','var');

if lye+lxe+uye+uxe==0 | isempty(lx) & isempty(ux) & isempty(ly) & isempty(uy)
    return
end

if uye==0 | isempty(uy)
    uy=ly;
end

if uxe==0 | isempty(ux)
    ux=lx;
end

for t=1:length(x)
    
if ~isempty(ux)
    %x errorbars
        PlotHandle(end+1)=line([x(t)-lx(t) x(t)+ux(t)],[y(t) y(t)],'color',errorcol);
        %(*) PlotHandle(end+1)=line([x(t)-lx(t) x(t)-lx(t)],[y(t)-yw y(t)+yw],'color',errorcol);  
        %(*) PlotHandle(end+1)=line([x(t)+ux(t) x(t)+ux(t)],[y(t)-yw y(t)+yw],'color',errorcol);
    end
    
    if ~isempty(uy)
    %y errorbars
        PlotHandle(end+1)=line([x(t) x(t)],[y(t)-ly(t) y(t)+uy(t)],'color',errorcol);
        %(*) PlotHandle(end+1)=line([x(t)-xw x(t)+xw],[y(t)-ly(t) y(t)-ly(t)],'color',errorcol);
        %(*) PlotHandle(end+1)=line([x(t)-xw x(t)+xw],[y(t)+uy(t) y(t)+uy(t)],'color',errorcol);
    end    
end

if ~KeepHold
    hold off
end
    