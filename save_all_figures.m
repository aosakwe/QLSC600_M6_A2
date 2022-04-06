function []=save_all_figures(outdir)
outdir=fullfile(outdir,'Figure');
if ~isfolder(outdir)
    mkdir(outdir)
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    savefig(FigHandle, fullfile(outdir, append(FigName,int2str(length(FigList) - iFig + 1), '.fig')));
    saveas(FigHandle, fullfile(outdir, append(FigName,int2str(length(FigList) - iFig + 1), '.png')));
    clf(FigHandle);
end
end