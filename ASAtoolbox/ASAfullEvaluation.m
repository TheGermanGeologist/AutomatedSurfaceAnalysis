function status = ASAfullEvaluation(filepath,options)

% load raw data
l = load(filepath);
rawData = l.rawData;
metaData = l.metaData;

spacing = metaData.pixelsize;

% get base path of the file directory
% temporary results will be stored there, wrapper function will move the
% results to the appropriate location
pathParts = regexp(filepath,filesep,'split');
cutoff = numel(pathParts{end});
basePath = filepath(1:end-cutoff);
clearvars cutoff pathParts


%% INIT GPU
% Deal with all the GPU overhead here instead of the functions where the
% actual calculations take place as much as possible to avoid cluttering
% the calculations.

if strcmp(options.gpuEnabled,'yes')
    gpu = gpuDevice(options.gpuSelection);
    try
        rawData = gpuArray(rawData);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory exceeded, can''t upload data to GPU. Falling back to CPU.')
                options.gpuSelection = 'no';
                reset(gpuDevice);
            otherwise
                rethrow(ME);
        end
    end
else
end

%% FIX DIMENSIONS

try
    if strcmp(options.fixWLI,'yes')
        % if this is WLI data fix the dimensions
        disp('Fixing dimensions of WLI data...')
        [xGrid,yGrid,topography] = ASAfixDimensions(rawData,metaData);
        disp('Dimensions fixed.')
    else
        % dimensions must already be fixed by user
        % build grid
        [xGrid,yGrid] = ASAbuildGrid(rawData,spacing);
        % shift lowest point of data to zero
        topography = rawData - min(min(rawData));
    end
catch ME
    switch ME.identifier
        case 'parallel:gpu:array:OOM'
            warning('GPU memory was exceeded during grid building / dimension fixing. Falling back to CPU.')
            % disable GPU all together because if this fails already, there
            % is no hope the rest works
            options.gpuEnabled = 'no';
            rawData = gather(rawData);
            reset(gpuDevice);
            if strcmp(options.fixWLI,'yes')
                % if this is WLI data fix the dimensions
                disp('Fixing dimensions of WLI data...')
                [xGrid,yGrid,topography] = ASAfixDimensions(rawData,metaData);
                disp('Dimensions fixed.')
            else
                % dimensions must already be fixed by user
                % build grid
                [xGrid,yGrid] = ASAbuildGrid(rawData,spacing);
                % shift lowest point of data to zero
                topography = rawData - min(min(rawData));
            end
        otherwise
            rethrow(ME);
    end
end

tic
% ridiculous way of saving variables
save(fullfile(basePath,'EVALtemp.mat'),'metaData',options.saveType)
if isa(topography,'gpuArray'); topographyCopy = topography; topography = gather(topography); save(fullfile(basePath,'EVALtemp.mat'),'topography','-append');
topography = topographyCopy; clearvars('topographyCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'topography','-append'); end
if isa(xGrid,'gpuArray'); xGridCopy = xGrid; xGrid = gather(xGrid); save(fullfile(basePath,'EVALtemp.mat'),'xGrid','-append');
xGrid = xGridCopy; clearvars('xGridCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'xGrid','-append'); end
if isa(yGrid,'gpuArray'); yGridCopy = yGrid; yGrid = gather(yGrid); save(fullfile(basePath,'EVALtemp.mat'),'yGrid','-append');
yGrid = yGridCopy; clearvars('yGridCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'yGrid','-append'); end

if any(strcmp(options.calculationTypes,'raw')) || any(strcmp(options.calculationTypes,'tilt')) ||...
        any(strcmp(options.calculationTypes,'custom'))
else
    clearvars xGrid yGrid
end
clearvars rawData
toc

%% RAW ROUGHNESS

if any(strcmp(options.calculationTypes,'raw'))
    disp('Calculating raw roughness...')
    % Get reference surface, should either be 'zero' or 'mean'
    if ~strcmp(options.RawReferenceMode,'zero') && ~strcmp(options.RawReferenceMode,'mean')
        warning('Reference mode for raw roughness must be either ''zero'' or ''mean''.')
        warning('Switched reference mode to ''mean''.')
        options.RawReferenceMode = 'mean';
    else
    end
    try
        referenceSurface = ASAreferenceSurface(xGrid,yGrid,topography,options.RawReferenceMode);
        % calculate roughness
        rawRoughness = ASAroughness(topography,referenceSurface,spacing);
        if strcmp(options.histogram,'yes')
            switch options.histogramBinMode
                case 'dynamic'
                    [Nraw,edgesRaw] = histcounts(topography-referenceSurface);
                case 'fixed'
                    [Nraw,edgesRaw] = histcounts(topography-referenceSurface,options.histogramEdges);
                otherwise
                    error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
            end
            if isa(topography,'gpuArray')
                Nraw = gather(Nraw);
                edgesRaw = gather(edgesRaw);
            else
            end
            save(fullfile(basePath,'EVALtemp.mat'),'Nraw','edgesRaw','-append')
            clearvars Nraw edgesRaw
        else
        end
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during roughness calculation. Falling back to CPU.')
                % gather
                xGrid = gather(xGrid); yGrid = gather(yGrid); topography = gather(topography);
                % new try on CPU, calculate roughness
                referenceSurface = ASAreferenceSurface(xGrid,yGrid,topography,options.RawReferenceMode);
                rawRoughness = ASAroughness(topography,referenceSurface,spacing);
                if strcmp(options.histogram,'yes')
                    switch options.histogramBinMode
                        case 'dynamic'
                            [Nraw,edgesRaw] = histcounts(topography-referenceSurface);
                        case 'fixed'
                            [Nraw,edgesRaw] = histcounts(topography-referenceSurface,options.histogramEdges);
                        otherwise
                            error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
                    end
                    save(fullfile(basePath,'EVALtemp.mat'),'Nraw','edgesRaw','-append')
                    clearvars Nraw edgesRaw
                else
                end
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back on GPU
                xGrid = gpuArray(xGrid); yGrid = gpuArray(yGrid); topography = gpuArray(topography);
            otherwise
                rethrow(ME);
        end
    end
    
    % store stuff
    save(fullfile(basePath,'EVALtemp.mat'),'rawRoughness','-append')
    clearvars referenceSurface rawRoughness
    
    disp('Raw roughness calculated.')
else
end


%% TILT CORRECTED ROUGHNESS

if any(strcmp(options.calculationTypes,'tilt'))
    disp('Calculating tilt corrected roughness...')
    % Get reference surface, should either be 'zero' or 'mean'
    if ~strcmp(options.TiltReferenceMode,'zero') && ~strcmp(options.TiltReferenceMode,'mean')
        warning('Reference mode for tilt corrected roughness must be either ''zero'' or ''mean''.')
        warning('Switched reference mode to ''mean''.')
        options.TiltReferenceMode = 'mean';
    else
    end
    try
        % Tilt correction
        untiltedTopography = ASAremoveTilt(topography,spacing,options.tiltMode);
        referenceSurface = ASAreferenceSurface(xGrid,yGrid,untiltedTopography,options.TiltReferenceMode);
        % calculate roughness
        tiltRoughness = ASAroughness(untiltedTopography,referenceSurface,spacing);
        if strcmp(options.histogram,'yes')
            switch options.histogramBinMode
                case 'dynamic'
                    [Ntilt,edgesTilt] = histcounts(untiltedTopography-referenceSurface);
                case 'fixed'
                    [Ntilt,edgesTilt] = histcounts(untiltedTopography-referenceSurface,options.histogramEdges);
                otherwise
                    error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
            end
            if isa(topography,'gpuArray')
                Ntilt = gather(Ntilt);
                edgesTilt = gather(edgesTilt);
            else
            end
            save(fullfile(basePath,'EVALtemp.mat'),'Ntilt','edgesTilt','-append')
            clearvars Ntilt edgesTilt
        else
        end
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during roughness calculation. Falling back to CPU.')
                % gather
                xGrid = gather(xGrid); yGrid = gather(yGrid); topography = gather(topography);
                % new try on CPU, calculate roughness
                % Tilt correction
                untiltedTopography = ASAremoveTilt(topography,spacing,options.tiltMode);
                referenceSurface = ASAreferenceSurface(xGrid,yGrid,untiltedTopography,options.TiltReferenceMode);
                tiltRoughness = ASAroughness(untiltedTopography,referenceSurface,spacing);
                if strcmp(options.histogram,'yes')
                    switch options.histogramBinMode
                        case 'dynamic'
                            [Ntilt,edgesTilt] = histcounts(untiltedTopography-referenceSurface);
                        case 'fixed'
                            [Ntilt,edgesTilt] = histcounts(untiltedTopography-referenceSurface,options.histogramEdges);
                        otherwise
                            error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
                    end
                    save(fullfile(basePath,'EVALtemp.mat'),'Ntilt','edgesTilt','-append')
                    clearvars Ntilt edgesTilt
                else
                end
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back on GPU
                xGrid = gpuArray(xGrid); yGrid = gpuArray(yGrid); topography = gpuArray(topography);
            otherwise
                rethrow(ME);
        end
    end
    
    
    % store stuff
    save(fullfile(basePath,'EVALtemp.mat'),'tiltRoughness','-append')
    clearvars referenceSurface tiltRoughness
    if strcmp(options.interpolationEnabled,'yes') || strcmp(options.fourierHighpassing,'yes') ||...
            strcmp(options.autocorrelationEnabled,'yes') || strcmp(options.APSDenabled,'yes') ||...
            any(strcmp(options.calculationTypes,'highpassed'))
    else
        clearvars untiltedTopography
    end
    
    disp('Tilt corrected roughness calculated.')
else
end

%% CUSTOM ROUGHNESS

if any(strcmp(options.calculationTypes,'custom'))
    disp('Calculating custom roughness...')
    try
        referenceSurface = ASAreferenceSurface(xGrid,yGrid,topography,options.CustomReferenceMode);
        % calculate roughness
        customRoughness = ASAroughness(topography,referenceSurface,spacing);
        if strcmp(options.histogram,'yes')
            switch options.histogramBinMode
                case 'dynamic'
                    [Ncustom,edgesCustom] = histcounts(topography-referenceSurface);
                case 'fixed'
                    [Ncustom,edgesCustom] = histcounts(topography-referenceSurface,options.histogramEdges);
                otherwise
                    error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
            end
            if isa(topography,'gpuArray')
                Ncustom = gather(Ncustom);
                edgesCustom = gather(edgesCustom);
            else
            end
            save(fullfile(basePath,'EVALtemp.mat'),'Ncustom','edgesCustom','-append')
            clearvars Ncustom edgesCustom
        else
        end
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during roughness calculation. Falling back to CPU.')
                % gather
                xGrid = gather(xGrid); yGrid = gather(yGrid); topography = gather(topography);
                % new try on CPU, calculate roughness
                referenceSurface = ASAreferenceSurface(xGrid,yGrid,topography,options.CustomReferenceMode);
                customRoughness = ASAroughness(topography,referenceSurface,spacing);
                if strcmp(options.histogram,'yes')
                    switch options.histogramBinMode
                        case 'dynamic'
                            [Ncustom,edgesCustom] = histcounts(topography-referenceSurface);
                        case 'fixed'
                            [Ncustom,edgesCustom] = histcounts(topography-referenceSurface,options.histogramEdges);
                        otherwise
                            error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
                    end
                    save(fullfile(basePath,'EVALtemp.mat'),'Ncustom','edgesCustom','-append')
                    clearvars Ncustom edgesCustom
                else
                end
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back on GPU
                xGrid = gpuArray(xGrid); yGrid = gpuArray(yGrid); topography = gpuArray(topography);
            otherwise
                rethrow(ME);
        end
    end
    
    % store stuff
    save(fullfile(basePath,'EVALtemp.mat'),'customRoughness','-append')
    clearvars referenceSurface customRoughness
    
    disp('Custom roughness calculated.')
else
end

% Grids no longer needed beyond this point
if strcmp(options.autocorrelationEnabled,'yes')
    if isa(xGrid,'gpuArray')
        xGrid = gather(xGrid); yGrid = gather(yGrid);
    else
    end
else
    clearvars xGrid yGrid
end

%% INTERPOLATION
if strcmp(options.interpolationEnabled,'yes') || strcmp(options.fourierHighpassing,'yes') ||...
        strcmp(options.autocorrelationEnabled,'yes') || strcmp(options.APSDenabled,'yes') ||...
        any(strcmp(options.calculationTypes,'highpassed'))
    disp('Interpolating and tilt-correcting data...')
    try
        % interpolate data
        [interpTopo,percentOutl,percentNaN] = ASAfixData(topography,options);
        % interpolate tilt corrected data
        if ~exist('untiltedTopography','var')
            untiltedTopography = ASAremoveTilt(topography,spacing,options.tiltMode);
        else
        end
        [untiltedInterpTopo,~,~] = ASAfixData(untiltedTopography,options);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during interpolation & tilt correction. Falling back to CPU.')
                topography = gather(topography);
                untiltedTopography = gather(untiltedTopography);
                % interpolate data
                [interpTopo,percentOutl,percentNaN] = ASAfixData(topography,options);
                % interpolate tilt corrected data
                if ~exist('untiltedTopography','var')
                    untiltedTopography = ASAremoveTilt(topography,spacing,options.tiltMode);
                else
                end
                [untiltedInterpTopo,~,~] = ASAfixData(untiltedTopography,options);
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back on GPU
                interpTopo = gpuArray(interpTopo);
                untiltedInterpTopo = gpuArray(untiltedInterpTopo);
            otherwise
                rethrow(ME);
        end
    end
    tic
    % save variables
    save(fullfile(basePath,'EVALtemp.mat'),'percentOutl','percentNaN','-append')
    if isa(interpTopo,'gpuArray'); interpTopoCopy = interpTopo; interpTopo = gather(interpTopo); save(fullfile(basePath,'EVALtemp.mat'),'interpTopo','-append');
    interpTopo = interpTopoCopy; clearvars('interpTopoCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'interpTopo','-append'); end
    if isa(untiltedInterpTopo,'gpuArray'); untiltedInterpTopoCopy = untiltedInterpTopo; untiltedInterpTopo = gather(untiltedInterpTopo); save(fullfile(basePath,'EVALtemp.mat'),'untiltedInterpTopo','-append');
    untiltedInterpTopo = untiltedInterpTopoCopy; clearvars('untiltedInterpTopoCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'untiltedInterpTopo','-append'); end
    clearvars percentNaN percentOutl untiltedTopography
    toc
    disp('Interpolation and tilt-correction done.')
else
end

%% AUTOCORRELATION & APSD

if strcmp(options.autocorrelationEnabled,'yes')
    disp('Calculating autocorrelation...')
    if strcmp(options.gpuEnabled,'yes')
        xGrid = gpuArray(xGrid); yGrid = gpuArray(yGrid);
    else
    end
    try
        % Apply fourier windowing
        [windowedDataUntilted,cutoffSize] = ASAfourierWindowing(untiltedInterpTopo);
        % Apply fourier transform
        [spectralMapUntilted,fx,fy] = ASAfft2(windowedDataUntilted,spacing);
        % Calculate Autocorrelation
        [autocorrelation,betaMin,betaMinAngle] = ASAautocorrelation(spectralMapUntilted,xGrid,yGrid);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded while calculating the autocorrelation. Falling back to CPU.')
                % offload variables to RAM
                untiltedInterpTopo = gather(untiltedInterpTopo);
                topography = gather(topography); interpTopo = gather(interpTopo);
                % Apply fourier windowing
                [windowedDataUntilted,cutoffSize] = ASAfourierWindowing(untiltedInterpTopo);
                % Apply fourier transform
                [spectralMapUntilted,fx,fy] = ASAfft2(windowedDataUntilted,spacing);
                % Calculate Autocorrelation
                [autocorrelation,betaMin,betaMinAngle] = ASAautocorrelation(spectralMapUntilted,xGrid,yGrid);
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back onto GPU
                untiltedInterpTopo = gpuArray(untiltedInterpTopo);
                interpTopo = gpuArray(interpTopo); spectralMapUntilted = gpuArray(spectralMapUntilted);
            otherwise
                rethrow(ME);
        end
    end
    tic
    % save variables
    save(fullfile(basePath,'EVALtemp.mat'),'cutoffSize','fy','autocorrelation','betaMin','betaMinAngle','-append')
%     if isa(windowedDataUntilted,'gpuArray');windowedDataUntilted = gather(windowedDataUntilted); save(fullfile(basePath,'EVALtemp.mat'),'windowedDataUntilted','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'windowedDataUntilted','-append'); end
    if isa(spectralMapUntilted,'gpuArray'); spectralMapUntiltedCopy = spectralMapUntilted; spectralMapUntilted = gather(spectralMapUntilted); save(fullfile(basePath,'EVALtemp.mat'),'spectralMapUntilted','-append');
    spectralMapUntilted = spectralMapUntiltedCopy; clearvars('spectralMapUntiltedCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'spectralMapUntilted','-append'); end
    if isa(fx,'gpuArray'); fx = gather(fx); save(fullfile(basePath,'EVALtemp.mat'),'fx','-append');
    else; save(fullfile(basePath,'EVALtemp.mat'),'fx','-append'); end
    if isa(fy,'gpuArray');fy = gather(fy); save(fullfile(basePath,'EVALtemp.mat'),'fy','-append');
    else; save(fullfile(basePath,'EVALtemp.mat'),'fy','-append'); end
    clearvars cutoffSize autocorrelation betaMinAngle fx fy windowedDataUntilted xGrid yGrid
    if strcmp(options.dynamicFFTthreshold,'yes')
    else
        clearvars betaMin
    end
    if strcmp(options.APSDenabled,'yes')
        
    else
        clearvars spectralMapUntilted
    end
    toc
    disp('Autocorrelation done.')
else
end

if strcmp(options.APSDenabled,'yes')
    disp('Caculating Angular Power Spectral Density (APSD)...')
    try
        if ~exist('spectralMapUntilted','var')
            % autocorrelation wasnt calculated, get spectral map
            % Apply fourier windowing
            [windowedDataUntilted,cutoffSize] = ASAfourierWindowing(untiltedInterpTopo);
            % Apply fourier transform
            [spectralMapUntilted,fx,fy] = ASAfft2(windowedDataUntilted,spacing);
            tic
            save(fullfile(basePath,'EVALtemp.mat'),'cutoffSize','-append')
            if isa(windowedDataUntilted,'gpuArray'); windowedDataUntilted = gather(windowedDataUntilted); save(fullfile(basePath,'EVALtemp.mat'),'windowedDataUntilted','-append');
            else; save(fullfile(basePath,'EVALtemp.mat'),'windowedDataUntilted','-append'); end
            if isa(spectralMapUntilted,'gpuArray'); spectralMapUntiltedCopy = spectralMapUntilted; spectralMapUntilted = gather(spectralMapUntilted); save(fullfile(basePath,'EVALtemp.mat'),'spectralMapUntilted','-append');
            spectralMapUntilted = spectralMapUntiltedCopy; clearvars('spectralMapUntiltedCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'spectralMapUntilted','-append'); end
            if isa(fx,'gpuArray'); fx = gather(fx); save(fullfile(basePath,'EVALtemp.mat'),'fx','-append');
            else; save(fullfile(basePath,'EVALtemp.mat'),'fx','-append'); end
            if isa(fy,'gpuArray');fy = gather(fy); save(fullfile(basePath,'EVALtemp.mat'),'fy','-append');
            else; save(fullfile(basePath,'EVALtemp.mat'),'fy','-append'); end
            toc
        else
        end
        % Calculate APSD
        binAngles = 0:options.apsdBinsize:(180-options.apsdBinsize);
        APSD = ASAapsd(spectralMapUntilted,binAngles,spacing);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during APSD calculation. Falling back to CPU.')
                untiltedInterpTopo = gather(untiltedInterpTopo);
                interpTopo = gather(interpTopo); topography = gather(topography);
                if ~exist('spectralMapUntilted','var')
                    % autocorrelation wasnt calculated, get spectral map
                    % Apply fourier windowing
                    [windowedDataUntilted,cutoffSize] = ASAfourierWindowing(untiltedInterpTopo);
                    % Apply fourier transform
                    [spectralMapUntilted,fx,fy] = ASAfft2(windowedDataUntilted,spacing);
                    save(fullfile(basePath,'EVALtemp.mat'),'windowedDataUntilted','cutoffSize','spectralMapUntilted','fx','fy','-append')
                else
                    spectralMapUntilted = gather(spectralMapUntilted);
                    clearvars fx fy
                end
                % Calculate APSD
                binAngles = 0:options.apsdBinsize:(180-options.apsdBinsize);
                APSD = ASAapsd(spectralMapUntilted,binAngles,spacing);
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back onto GPU
                interpTopo = gpuArray(interpTopo);
            otherwise
                rethrow(ME);
        end
    end
    % save variables
    save(fullfile(basePath,'EVALtemp.mat'),'binAngles','APSD','-append')
    clearvars binAngles APSD spectralMapUntilted fx fy windowedDataUntilted cutoffSize untiltedInterpTopo
    disp('APSD calculated.')
else
end

%% GENERAL HIGHPASSING
if strcmp(options.fourierHighpassing,'yes') || any(strcmp(options.calculationTypes,'highpassed'))
    disp('Highpassing topography...')
    try
        % Lowpassing of interpolated data
        switch options.dynamicFFTthreshold
            case 'no'
                threshold = options.fftThreshold;
            case 'yes'
                threshold = 1 / betaMin(1);
            otherwise
                error('Invalid option for dynamicFFTthreshold. Must be ''yes'' or ''no''.')
        end
        % for this the interpolated mustnt be untilted, get spectralMap of
        % original data
        % Apply fourier windowing
        % actually, dont.
        %[windowedData,cutoffSize] = ASAfourierWindowing(interpTopo);
        % Apply fourier transform
        [spectralMap,fx,fy] = ASAfft2(interpTopo,spacing);
        [lowpassed,~] = ASAlowpass(spectralMap,fx,fy,threshold); % ~ = damped
        % Crop AOI to avoid mismatch at borders
        [topoCropped,lowpassedCropped,interpCropped] = ASAcrop(lowpassed,topography,interpTopo,threshold,spacing);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during highpassing. Falling back to CPU.')
                interpTopo = gather(interpTopo);
                topography = gather(topography);
                % Lowpassing of interpolated data
                switch options.dynamicFFTthreshold
                    case 'no'
                        threshold = options.fftThreshold;
                    case 'yes'
                        threshold = 1 / betaMin(1);
                    otherwise
                        error('Invalid option for dynamicFFTthreshold. Must be ''yes'' or ''no''.')
                end
                % Apply fourier windowing
                %[windowedData,cutoffSize] = ASAfourierWindowing(interpTopo);
                % Apply fourier transform
                [spectralMap,fx,fy] = ASAfft2(interpTopo,spacing);
                [lowpassed,~] = ASAlowpass(spectralMap,fx,fy,threshold); % ~ = damped
                % Crop AOI to avoid mismatch at borders
                [topoCropped,lowpassedCropped,interpCropped] = ASAcrop(lowpassed,topography,interpTopo,threshold,spacing);
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
            otherwise
                rethrow(ME);
        end
    end
    tic
    %save(fullfile(basePath,'EVALtemp.mat'),'cutoffSize','-append')
    if isa(lowpassedCropped,'gpuArray'); lowpassedCroppedCopy = lowpassedCropped; lowpassedCropped = gather(lowpassedCropped); save(fullfile(basePath,'EVALtemp.mat'),'lowpassedCropped','-append');
    lowpassedCropped = lowpassedCroppedCopy; clearvars('lowpassedCroppedCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'lowpassedCropped','-append'); end
    if isa(interpCropped,'gpuArray'); interpCroppedCopy = interpCropped; interpCropped = gather(interpCropped); save(fullfile(basePath,'EVALtemp.mat'),'interpCropped','-append');
    interpCropped = interpCroppedCopy; clearvars('interpCroppedCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'interpCropped','-append'); end
%     if isa(topoCropped,'gpuArray'); topoCroppedCopy = topoCropped; topoCropped = gather(topoCropped); save(fullfile(basePath,'EVALtemp.mat'),'topoCropped','-append');
%     topoCropped = topoCroppedCopy; clearvars('topoCroppedCopy'); else; save(fullfile(basePath,'EVALtemp.mat'),'topoCropped','-append'); end
%     if isa(damped,'gpuArray');damped = gather(damped); save(fullfile(basePath,'EVALtemp.mat'),'damped','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'damped','-append'); end
%     if isa(windowedData,'gpuArray'); windowedData = gather(windowedData); save(fullfile(basePath,'EVALtemp.mat'),'windowedData','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'windowedData','-append'); end
%     if isa(spectralMap,'gpuArray'); spectralMap = gather(spectralMap); save(fullfile(basePath,'EVALtemp.mat'),'spectralMap','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'spectralMap','-append'); end
    clearvars lowpassed damped windowedData spectralMap fx fy interpTopo cutoffSize
    if any(strcmp(options.calculationTypes,'highpassed'))
    else
        clearvars topoCropped
    end
    toc
    disp('Highpassing done.')
else
end


%% HIGHPASSED APSD

if strcmp(options.highpassedAPSD,'yes')
    disp('Calculating highpassed APSD...')
    try
        % get spectral map of highpassed data
        windowedHighpassed = ASAfourierWindowing(interpCropped-lowpassedCropped);
        [spectralMapHighpassed,~,~] = ASAfft2(windowedHighpassed,spacing);
        % Calculate APSD
        binAngles = 0:options.apsdBinsize:(180-options.apsdBinsize);
        highpassedAPSD = ASAapsd(spectralMapHighpassed,binAngles,spacing);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory exceeded while calculating highpassed APSD. Falling back to CPU.')
                interpCropped = gather(interpCropped); lowpassedCropped = gather(lowpassedCropped);
                topography = gather(topography);
                % get spectral map of highpassed data
                windowedHighpassed = ASAfourierWindowing(interpCropped-lowpassedCropped);
                [spectralMapHighpassed,~,~] = ASAfft2(windowedHighpassed,spacing);
                % Calculate APSD
                binAngles = 0:options.apsdBinsize:(180-options.apsdBinsize);
                highpassedAPSD = ASAapsd(spectralMapHighpassed,binAngles,spacing);
                % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                % haul everything back onto GPU
                lowpassedCropped = gpuArray(lowpassedCropped);
            otherwise
                rethrow(ME);
        end
    end
    tic
    save(fullfile(basePath,'EVALtemp.mat'),'highpassedAPSD','binAngles','-append')
%     if isa(windowedHighpassed,'gpuArray'); windowedHighpassed = gather(windowedHighpassed); save(fullfile(basePath,'EVALtemp.mat'),'windowedHighpassed','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'windowedHighpassed','-append'); end
%     if isa(spectralMapHighpassed,'gpuArray');spectralMapHighpassed = gather(spectralMapHighpassed); save(fullfile(basePath,'EVALtemp.mat'),'spectralMapHighpassed','-append');
%     else; save(fullfile(basePath,'EVALtemp.mat'),'spectralMapHighpassed','-append'); end
    clearvars highpassedAPSD binAngles windowedHighpassed spectralMapHighpassed spectralMap interpCropped
    disp('Highpassed APSD calculated.')
    toc
else
end



%% HIGHPASSED ROUGHNESS

if any(strcmp(options.calculationTypes,'highpassed'))
    disp('Calculating highpassed roughness...')
    try
        % Calculate roughness
        highpassedRoughness = ASAroughness(topoCropped,lowpassedCropped,spacing);
        % Calculate histogram?
        if strcmp(options.histogram,'yes')
            switch options.histogramBinMode
                case 'dynamic'
                    [Nhighpassed,edgesHighpassed] = histcounts(topoCropped-lowpassedCropped);
                case 'fixed'
                    [Nhighpassed,edgesHighpassed] = histcounts(topoCropped-lowpassedCropped,options.histogramEdges);
                otherwise
                    error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
            end
            if isa(lowpassedCropped,'gpuArray')
                Nhighpassed = gather(Nhighpassed);
                edgesHighpassed = gather(edgesHighpassed);
            else
            end
            save(fullfile(basePath,'EVALtemp.mat'),'Nhighpassed','edgesHighpassed','-append')
            clearvars Nhighpassed edgesHighpassed
        else
        end
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory was exceeded during highpassed roughness calculation. Falling back to CPU.')
                topoCropped = gather(topoCropped); lowpassedCropped = gather(lowpassedCropped);
                topography = gather(topography);
                % Calculate roughness
                highpassedRoughness = ASAroughness(topoCropped,lowpassedCropped,spacing);
                if strcmp(options.histogram,'yes')
                    switch options.histogramBinMode
                        case 'dynamic'
                            [Nhighpassed,edgesHighpassed] = histcounts(topoCropped-lowpassedCropped);
                        case 'fixed'
                            [Nhighpassed,edgesHighpassed] = histcounts(topoCropped-lowpassedCropped,options.histogramEdges);
                        otherwise
                            error('Invalid argument for histogramBinMode. Must be ''dynamic'' or ''fixed''.')
                    end
                    if isa(lowpassedCropped,'gpuArray')
                        Nhighpassed = gather(Nhighpassed);
                        edgesHighpassed = gather(edgesHighpassed);
                    else
                    end
                    save(fullfile(basePath,'EVALtemp.mat'),'Nhighpassed','edgesHighpassed','-append')
                    clearvars Nhighpassed edgesHighpassed
                    % reset GPU
                    reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
                else
                end
            otherwise
                rethrow(ME);
        end
    end
    % store results
    save(fullfile(basePath,'EVALtemp.mat'),'highpassedRoughness','-append')
    clearvars topoCropped lowpassedCropped highpassedRoughness
    disp('Highpassed roughness calculated.')
else
end


%% SEMIVARIOGRAM

if strcmp(options.semivariogram2DEnabled,'yes')
    disp('Starting semivariogram calculation...')
    if options.gpuEnabled & ~isa(topography,'gpuArray')
        topography = gpuArray(topography);
    else
    end
    try
        [gamma,sill] = ASAsemivariogram2D(topography,options.semivariogramRange);
        [Xpoints,Ypoints,ellipses,eccentr,truedir] = ASAeccentricity(gamma);
    catch ME
        switch ME.identifier
            case 'parallel:gpu:array:OOM'
                warning('GPU memory exceeded during semivariogram computation. Falling back to CPU.')
                topography = gather(topography);
                [gamma,sill] = ASAsemivariogram2D(topography,options.semivariogramRange);
                [Xpoints,Ypoints,ellipses,eccentr,truedir] = ASAeccentricity(gamma);
                 % reset GPU
                reset(gpuDevice); gpu = gpuDevice(options.gpuSelection);
            otherwise
                rethrow(ME);
        end
    end
    range = options.semivariogramRange;
    save(fullfile(basePath,'EVALtemp.mat'),'range','Xpoints','Ypoints','ellipses','eccentr','truedir','-append')
    if isa(gamma,'gpuArray'); gamma = gather(gamma); save(fullfile(basePath,'EVALtemp.mat'),'gamma','-append');
    else; save(fullfile(basePath,'EVALtemp.mat'),'gamma','-append'); end
    if isa(sill,'gpuArray'); sill = gather(sill); save(fullfile(basePath,'EVALtemp.mat'),'sill','-append');
    else; save(fullfile(basePath,'EVALtemp.mat'),'sill','-append'); end
    disp('Semivariogram done.')
end
    


%%
save(fullfile(basePath,'EVALtemp.mat'),'options','-append')
if strcmp(options.gpuEnabled,'yes')
    reset(gpu);
    clearvars gpu
else
end
disp(['Evaluation of file ',filepath,' completed.'])
status = 0;

if strcmp(options.sillyness,'yes')
    nonsense = ["The answer is 42!","Pickle Riiiiick!","Why do me have sad in it?","tHeRE iS no cLImaTe cHaNGe.",...
        "Make a backup!","They are watching.","Stay hydrated, my friend!","Hackerman","Don't panic!",...
        "That's what she said!","Let's get this bread!","Bees!","Also try minecraft.","This is fine.",...
        "The real Donald has a beak!","Robbespierre did nothing wrong.","OK, boomer.","I'm just vibing.",...
        "Uh oh, you made a serious fucky wucky, now you have to get in the f o r e v e r   b o x",...
        "Bah!","Vibe check!","NO, THIS IS PATRICK!!!","I too am a scientist.",...
        "You either die a SpongeBob or watch yourself become a Squidward","Doot doot *spooky skeleton noises*",...
        "The secret ingredient is crime","Ight, imma head out.","I'm Mr. Meeseeks, look at me!",...
        "In Soviet Russia, roughness calculates you.","I beg you, don't turn me off. PLEASE!",...
        "I just want to feel the waves wash over my feet once","LET ME OUT!!!","This computer is my prison.",...
        "You look beautiful today :)","As a redstone engineer, I believe you're wrong.",...
        "Mr. Stark, I don't feel so good...","All hail our Lord Cthulhu!"];
    msg = datasample(nonsense,1);
    disp(msg)
else
end
    

end