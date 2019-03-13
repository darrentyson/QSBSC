cd ~/Dropbox/cancer-rna-fish/ 
% change this to the location of the repo on your computer

T = readtable('dentistData/WM9_noDrug_20150618.txt');

Xpos = T.Xpos;
Ypos = T.Ypos;
AXL = T.AXL;
cellID = T.cellID;


colormap default
cmap = colormap('jet');

cmap = cmap(10:55,:);


scatter(Xpos,Ypos,100,AXL,'.');
colormap(cmap);
axis equal
%set(gcf,'PaperPositionMode','manual');
print -depsc2 graphs/exampleDentistPlots/WM9_noDrug_20150618_dentistDataPlotAXL.eps

% 257 spots in Cell ID 1410

idx = find(cellID == 1410);

hold on;
plot(Xpos(idx),Ypos(idx),'ro','markersize',15);
print -depsc2 graphs/exampleDentistPlots/WM9_noDrug_20150618_dentistDataPlotAXL_CellCallout.eps

colorbar
print -depsc2 graphs/exampleDentistPlots/WM9_noDrug_20150618_dentistDataPlotAXL_colorbar.eps



