function ASAgenReport(reportpath,resultsStruct,evalname)

ln1 = ['Report for ASA evaluation of scan ',evalname,' :'];
lnsep = '________________________________________________';
ln3 = ['Evaluation performed on ',datestr(datetime)];
% lineseparator
ln4 = ['>>> METADATA<<<'];
ln5 = ['Scan type: ',];
ln6 = ['Material type: ',];
ln7 = ['Deformation type: ',];
ln8 = ['Scan dimensions (x,y): ',num2str(fliplr(size(resultsStruct.topography)))];
% lineseparator
ln9 = '>>>OPTIONS<<<';
ln10 = 




end