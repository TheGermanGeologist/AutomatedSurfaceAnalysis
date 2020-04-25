function ASAfigures(evalfilepath,rawfilepath,savefolderpath,evalname,options)

close all

l = load(evalfilepath);
r = load(rawfilepath);
figurepath = fullfile(savefolderpath,'figures');
mkdir(figurepath)

%% RAW topography
if isfield(l,'topography')
    v = figure('Name','Raw data','Position',[0 0 720 720]);
    f = imagesc(l.xGrid(1,:),l.yGrid(:,1),l.topography);
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    %colormap(demcmapc)
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    set(f,'AlphaData',~isnan(l.topography))
    title([evalname,': Raw data'])
    xlim([min(l.xGrid(1,:)) max(l.xGrid(1,:))])
    ylim([min(l.yGrid(:,1)) max(l.yGrid(:,1))])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'rawdata'),options.imageFormats(i))
        end
    else
    end
else
end

%% INTERPOLATED topography
if isfield(l,'interpTopo')
    v =  figure('Name','Interpolated data','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,:),l.yGrid(:,1),l.interpTopo)
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    %colormap(demcmapc)
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    title([evalname,': Interpolated data'])
    xlim([min(l.xGrid(1,:)) max(l.xGrid(1,:))])
    ylim([min(l.yGrid(:,1)) max(l.yGrid(:,1))])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'interpolated'),options.imageFormats(i))
        end
    else
    end
else
end

%% UNTILTED INTERPOLATED topography
if isfield(l,'untiltedInterpTopo')
    v = figure('Name','Interpolated & tilt removed','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,:),l.yGrid(:,1),l.untiltedInterpTopo)
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    %colormap(demcmapc)
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    title([evalname,': Interpolated & tilt corrected data'])
    xlim([min(l.xGrid(1,:)) max(l.xGrid(1,:))])
    ylim([min(l.yGrid(:,1)) max(l.yGrid(:,1))])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'interpolated-untilted'),options.imageFormats(i))
        end
    else
    end
else
end

%% PROFILE LINES

if isfield(l,'interpTopo')
    [ySize,xSize] = size(l.interpTopo);
    
    yAxisVec = l.yGrid(:,1);
    xAxisVec = l.xGrid(1,:);
    
    v = figure('Name','X profile lines','Position',[0 0 720 360]);
    plot(xAxisVec,l.interpTopo(round(ySize*0.25),:))
    hold on
    plot(xAxisVec,l.interpTopo(round(ySize*0.5),:))
    plot(xAxisVec,l.interpTopo(round(ySize*0.75),:))
    xlabel('x distance [mm]')
    ylabel('Elevation [µm]')
    title([evalname,': Profile lines along the X-direction'])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'xProfiles'),options.imageFormats(i))
        end
    else
    end
    
    v = figure('Name','Y profile lines','Position',[0 0 720 360]);
    plot(yAxisVec,l.interpTopo(:,round(xSize*0.25)))
    hold on
    plot(yAxisVec,l.interpTopo(:,round(xSize*0.5)))
    plot(yAxisVec,l.interpTopo(:,round(xSize*0.75)))
    xlabel('y distance [mm]')
    ylabel('Elevation [µm]')
    title([evalname,': Profile lines along the Y-direction'])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'yProfiles'),options.imageFormats(i))
        end
    else
    end
    
    v = figure('Name','Profile lines map','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,:),l.yGrid(:,1),l.interpTopo)
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    title([evalname,': Interpolated data - profile lines'])
    xlim([min(l.xGrid(1,:)) max(l.xGrid(1,:))])
    ylim([min(l.yGrid(:,1)) max(l.yGrid(:,1))])
    hold on;
    line([l.xGrid(1,1),l.xGrid(1,xSize)],[l.yGrid(round(ySize*0.25),1),l.yGrid(round(ySize*0.25),1)], 'Color', 'k');
    line([l.xGrid(1,1),l.xGrid(1,xSize)],[l.yGrid(round(ySize*0.5),1),l.yGrid(round(ySize*0.5),1)], 'Color', 'k');
    line([l.xGrid(1,1),l.xGrid(1,xSize)],[l.yGrid(round(ySize*0.75),1),l.yGrid(round(ySize*0.75),1)], 'Color', 'k');
    line([l.xGrid(1,round(xSize*0.25)),l.xGrid(1,round(xSize*0.25))],[l.yGrid(1,1),l.yGrid(ySize,1)], 'Color', 'k');
    line([l.xGrid(1,round(xSize*0.5)),l.xGrid(1,round(xSize*0.5))],[l.yGrid(1,1),l.yGrid(ySize,1)], 'Color', 'k');
    line([l.xGrid(1,round(xSize*0.75)),l.xGrid(1,round(xSize*0.75))],[l.yGrid(1,1),l.yGrid(ySize,1)], 'Color', 'k');
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'interpolated-profiles'),options.imageFormats(i))
        end
    else
    end
else
end

%% FOURIER DECOMPOSITION
if isfield(l,'lowpassedCropped')
    [croppedY,croppedX] = size(l.lowpassedCropped);
    [originalY,originalX] = size(l.interpTopo);
    paddingX = round((originalX-croppedX)/2);
    paddingY = round((originalY-croppedY)/2);
    
    v = figure('Name','Low-passed topography','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,paddingX:(end-paddingX)),l.yGrid(paddingY:(end-paddingY),1),l.lowpassedCropped)
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    %colormap(demcmapc)
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    title([evalname,': Low-frequency topography component'])
    xlim([min(l.xGrid(1,paddingX:(end-paddingX))) max(l.xGrid(1,paddingX:(end-paddingX)))])
    ylim([min(l.yGrid(paddingY:(end-paddingY),1)) max(l.yGrid(paddingY:(end-paddingY),1))])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'lowpassed'),options.imageFormats(i))
        end
    else
    end
    
    v = figure('Name','High-passed topography','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,paddingX:(end-paddingX)),l.yGrid(paddingY:(end-paddingY),1),l.interpCropped-l.lowpassedCropped)
    axis equal
    xlabel('x direction [mm]')
    ylabel('y direction [mm]')
    %colormap(demcmapc)
    c = colorbar;
    ylabel(c,'Elevation [µm]')
    title([evalname,': High-frequency topography component'])
    xlim([min(l.xGrid(1,paddingX:(end-paddingX))) max(l.xGrid(1,paddingX:(end-paddingX)))])
    ylim([min(l.yGrid(paddingY:(end-paddingY),1)) max(l.yGrid(paddingY:(end-paddingY),1))])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'highpassed'),options.imageFormats(i))
        end
    else
    end
else
end

if isfield(l,'spectralMap')
    v = figure('Name','2D Fourier transformation','Position',[0 0 720 720]);
    imagesc(l.fx,l.fy,real(log(abs(l.spectralMapUntilted))))
    c = colorbar;
    axis equal
    xlabel('x-frequency [cycles/mm]')
    ylabel('y-frequency [cycles/mm]')
    ylabel(c,'Fourier coefficient [µm^2]')
    title([evalname,': 2D Fourier tranformation'])
    xlim([min(l.fx) max(l.fx)])
    ylim([min(l.fy) max(l.fy)])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'fourier'),options.imageFormats(i))
        end
    else
    end
    
    v = figure('Name','2D Fourier transformation cropped','Position',[0 0 720 720]);
    imagesc(l.fx,l.fy,real(log(abs(l.spectralMapUntilted))))
    c = colorbar;
    axis equal
    xlabel('x-frequency [cycles/mm]')
    ylabel('y-frequency [cycles/mm]')
    ylabel(c,'Fourier coefficient [µm^2]')
    title([evalname,': 2D Fourier tranformation'])
    xlim([-15 15])
    ylim([-15 15])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'fourier-zoomed'),options.imageFormats(i))
        end
    else
    end
else
end


%% AUTOCORRELATION

if isfield(l,'autocorrelation')
    v = figure('Name','Autocorrelation','Position',[0 0 720 720]);
    imagesc(l.xGrid(1,:)-max(l.xGrid(1,:))/2,l.yGrid(:,1)-max(l.yGrid(:,1))/2,real(abs(l.autocorrelation)))
    axis equal
    xlabel('shift in x direction [mm]')
    ylabel('shift in y direction [mm]')
    c = colorbar;
    colormap jet
    ylabel(c,'ACF')
    title([evalname,': Autocorrelation'])
    hold on
    [autoGridX,autoGridY] = meshgrid(l.xGrid(1,:)-max(l.xGrid(1,:))/2,l.yGrid(:,1)-max(l.yGrid(:,1))/2);
    x = autoGridX(l.autocorrelation < (max(max(l.autocorrelation)*0.1+max(max(l.autocorrelation)*0.005))) & l.autocorrelation > (max(max(l.autocorrelation)*0.1-max(max(l.autocorrelation)*0.005))));
    y = autoGridY(l.autocorrelation < (max(max(l.autocorrelation)*0.1+max(max(l.autocorrelation)*0.005))) & l.autocorrelation > (max(max(l.autocorrelation)*0.1-max(max(l.autocorrelation)*0.005))));
    plot(x,y,'.k')
    xlim([min(l.xGrid(1,:)-max(l.xGrid(1,:))/2) max(l.xGrid(1,:)-max(l.xGrid(1,:))/2)])
    ylim([min(l.yGrid(:,1)-max(l.yGrid(:,1))/2) max(l.yGrid(:,1)-max(l.yGrid(:,1))/2)])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'autocorrelation'),options.imageFormats(i))
        end
    else
    end
else
end

%% SEMIVARIOGRAM
if isfield(l,'gamma')
    rangeVec = -l.range:l.range;
    rangeScaled = rangeVec .* r.metaData.pixelsize .*1e3;
    %             [xGridvario,yGridvario] = meshgrid(rangeScaled,rangeScaled);
    
    v = figure('Name','Semivariogram');
    imagesc(rangeScaled,rangeScaled,l.gamma)
    xlabel('x-shift [µm]')
    ylabel('y-shift [µm]')
    title([evalName,': Semivariogram'])
    xlim([min(rangeScaled) max(rangeScaled)])
    ylim([min(rangeScaled) max(rangeScaled)])
    c = colorbar;
    % FIX UNIT OF COLORBAR!!!
    ylabel(c,'semivariance [µm^2]')
    axis equal
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'semivariogram'),options.imageFormats(i))
        end
    else
    end
    
    v = figure('Name','Semivariogram with isolines','Position',[0 0 720 720]);
    imagesc(-l.range:l.range,-l.range:l.range,l.gamma)
    hold on
    plot(l.Xpoints{1},l.Ypoints{1},'.k','LineWidth',1.7)
    plot(l.Xpoints{2},l.Ypoints{2},'.k','LineWidth',1.7)
    plot(l.Xpoints{3},l.Ypoints{3},'.k','LineWidth',1.7)
    xlabel('x-shift [pixel]')
    ylabel('y-shift [pixel]')
    title([evalname,': Semivariogram'])
    xlim([min(-l.range:l.range) max(-l.range:l.range)])
    ylim([min(-l.range:l.range) max(-l.range:l.range)])
    c = colorbar;
    % FIX UNIT OF COLORBAR!!!
    ylabel(c,'semivariance [µm^2]')
    axis equal
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'semivariogram-isolines'),options.imageFormats(i))
        end
    else
    end
else
end


%% APSD
if isfield(l,'APSD')
    apsdVec = [l.APSD l.APSD]';
    [m,n] = size(apsdVec);
    mm = 4*m;
    r = zeros(mm,n);
    r(2:4:mm,:) = apsdVec;
    r(3:4:mm,:) = apsdVec;
    zz = deg2rad((0-options.apsdBinsize/2):options.apsdBinsize:(360-options.apsdBinsize/2));
    t = zeros(mm,1);
    t(2:4:mm) = zz(1:m);
    t(3:4:mm) = zz(2:m+1);
    
    v = figure('Name','APSD polar plot','Position',[0 0 720 720]);
    h = polar(t,r);
    [a,b] = pol2cart(t,r);
    A = reshape(a,4,numel(apsdVec));
    B = reshape(b,4,numel(apsdVec));
    patch(A,B,[0.2 0.2 0.5])
    title([evalname,': Angular PSD'])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'apsd'),options.imageFormats(i))
        end
    else
    end
else
end


%% HIGHPASSED APSD
if isfield(l,'highpassedAPSD')
    apsdVec = [l.highpassedAPSD l.highpassedAPSD]';
    [m,n] = size(apsdVec);
    mm = 4*m;
    r = zeros(mm,n);
    r(2:4:mm,:) = apsdVec;
    r(3:4:mm,:) = apsdVec;
    zz = deg2rad((0-options.apsdBinsize/2):options.apsdBinsize:(360-options.apsdBinsize/2));
    t = zeros(mm,1);
    t(2:4:mm) = zz(1:m);
    t(3:4:mm) = zz(2:m+1);
    
    v = figure('Name','APSD polar plot - high-passed','Position',[0 0 720 720]);
    h = polar(t,r);
    [a,b] = pol2cart(t,r);
    A = reshape(a,4,numel(apsdVec));
    B = reshape(b,4,numel(apsdVec));
    patch(A,B,[0.2 0.2 0.5])
    title([evalname,': Angular PSD - Highpassed'])
    if strcmp(options.autoSaveImages,'yes')
        for i = 1:numel(options.imageFormats)
            saveas(v,fullfile(figurepath,'aspd-highpassed'),options.imageFormats(i))
        end
    else
    end
else
end

%%
close all


end