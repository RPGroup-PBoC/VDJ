function StandardFigure(PlotHandle,AxisHandle)

%Make it such that it can identify things like '-o'.

for i=1:length(PlotHandle)
    %if ~(get(PlotHandle(i),'LineStyle')=='none')
	
    %Line style
    LineStyle=get(PlotHandle(i),'LineStyle');
    if (~isempty(strmatch(LineStyle,'-')))|...
            (~isempty(strmatch(LineStyle,'--')))|...
            (~isempty(strmatch(LineStyle,'-.')))|...
            (~isempty(strmatch(LineStyle,'-o')))
        set(PlotHandle(i),'LineWidth',1)
    end
    

    

    %Marker style
    if get(PlotHandle(i),'Marker')=='.'
        set(PlotHandle(i),'MarkerSize',15)
        set(PlotHandle(i),'LineWidth',1)
    elseif (get(PlotHandle(i),'Marker')=='d')|...
        (strcmp(get(PlotHandle(i),'Marker'),'square'))|...
        (get(PlotHandle(i),'Marker')=='o')
        set(PlotHandle(i),'MarkerSize',5)
        set(PlotHandle(i),'LineWidth',1)
        if ~isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'none'))
            set(PlotHandle(i),'MarkerFaceColor','w')
        end
    end

    %Color section
    if isempty(strmatch(get(PlotHandle(i),'Color'),'auto'))
        ChangeColor(PlotHandle(i),'Color')
    end

    if isempty(strmatch(get(PlotHandle(i),'MarkerEdgeColor'),'auto'))
        ChangeColor(PlotHandle(i),'MarkerEdgeColor')
    end

    if isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'auto'))
        ChangeColor(PlotHandle(i),'MarkerFaceColor')
    end

%     
%     
%     for j=1:length(ColorProperty)
%         
%         
%             
%         
%         [i,j]
%         
%         
%         CurrentColor=get(PlotHandle(i),ColorProperty{j});
%         
%         if strmatch(CurrentColor,'auto')
%             CurrentColor=get(PlotHandle(i),'MarkerFaceColor');
%         end
%            
%       
%     end
            
            
            
            
    
    
    
%     %Point style
%     if isempty(strmatch(get(PlotHandle(i),'LineStyle'),'none'))
%         set(PlotHandle(i),'LineWidth',1)
%     end
    
    %Marker style
%     if get(PlotHandle(i),'Marker')=='.'
%         set(PlotHandle(i),'MarkerSize',15)
%         set(PlotHandle(i),'LineWidth',1)
% 	elseif (get(PlotHandle(i),'Marker')=='d')|...
%             (get(PlotHandle(i),'Marker')=='s')|...
%             (get(PlotHandle(i),'Marker')=='o')
%         set(PlotHandle(i),'MarkerSize',5)
%         set(PlotHandle(i),'LineWidth',1)
%     end
end


set(get(AxisHandle,'XLabel'),'FontSize',14,'FontName','Lucida Sans')
set(get(AxisHandle,'YLabel'),'FontSize',14,'FontName','Lucida Sans')
set(get(AxisHandle,'ZLabel'),'FontSize',14,'FontName','Lucida Sans')
set(AxisHandle,'FontSize',12,'FontName','Lucida Sans')
% set(AxisHandle,...
%     'Box','off',...
%     'Color',[0.8510,0.8510,0.8510],...
%     'XColor','w','YColor','w','ZColor','w',...
%     'TickLength',[0.02,0.05])
