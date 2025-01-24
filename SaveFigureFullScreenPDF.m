% Save figure to full-screen pdf

function SaveFigureFullScreenPDF(FigHandle,OutputFileName,size_fig)

%     set(FigHandle,'Units','Inches');
%     pos = get(FigHandle,'Position');
%     set(FigHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
%     
%     print(OutputFileName,'-dpdf','-r0');

%    h=gcf;
    if nargin == 2
        size = [20 10];
    else
        size = size_fig;
    end
    set(FigHandle,'PaperSize',size);
    print(FigHandle,OutputFileName,'-fillpage','-dpdf')
    
end