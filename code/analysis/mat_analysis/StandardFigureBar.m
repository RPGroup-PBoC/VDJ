%Steph 9/10

%Modification of Hernan's StandardFigure to work with barseries objects.

function StandardFigureBar(PlotHandle,AxisHandle)

for i=1:length(PlotHandle)

    %Color section
    if isempty(strmatch(get(PlotHandle(i),'FaceColor'),'auto'))
        ChangeColor(PlotHandle(i),'FaceColor')
    end
    
    if isempty(strmatch(get(PlotHandle(i),'EdgeColor'),'auto'))
        ChangeColor(PlotHandle(i),'EdgeColor')
    end


end

set(get(AxisHandle,'XLabel'),'FontSize',14,'FontName','Lucida Sans')
set(get(AxisHandle,'YLabel'),'FontSize',14,'FontName','Lucida Sans')
set(get(AxisHandle,'ZLabel'),'FontSize',14,'FontName','Lucida Sans')
set(AxisHandle,'FontSize',12,'FontName','Lucida Sans')