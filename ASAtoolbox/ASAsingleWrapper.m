function status = ASAsingleWrapper(paths,evalname,options)

status = -1;

%% Technical stuff

% Check for free disk space before doing anything else
savepathParts = regexp(paths.savepath,filesep,'split');
if ispc
    % Windows free disk space check
    [cmdstatus,volumeInfo] = system(['WMIC LOGICALDISK GET FreeSpace,Name,Size | find /i "',savepathParts{1},'"']);
    if cmdstatus ~= 0
        warning('Could not check free disk space. Continuing.')
    else
    end
    volumeInfoSplit = strsplit(volumeInfo);
    freeDiskSpace = str2double(volumeInfoSplit{1})/(1024^3);
    if freeDiskSpace < options.minimumDisk
        warning('Not enough disk space remaining. Exiting.')
        status = -2;
        return
    else
    end
elseif ismac
    warning('Check for free disk space not implemented on Mac OS. Continuing.')
elseif isunix && ~ismac
    warning('Check for free disk space not implemented on Linux. Continuing.')
else
    warning('Unknown OS - cannot check for free disk space. Continuing.')
end

% check for GPU if option selected
if strcmp(options.gpuEnabled,'yes')
    if gpuDeviceCount == 0
        warning('No compatible GPU found, falling back to CPU.')
        options.gpuEnabled = 'no';
    elseif gpuDeviceCount == 1
        gpu = gpuDevice;
        options.gpuSelection = 1;
        if gpu.DeviceSupported ~= 1
            warning('GPU not supported, falling back to CPU.')
            options.gpuEnabled = 'no';
        else
        end
        if gpu.TotalMemory/(1024^3) < 4
            warning('GPU has less than the recommended 4 GB of memory.')
        else
        end
    elseif gpuDeviceCount > 1
        disp('Multiple GPUs detected, trying to find best one.')
        gpuSupported = zeros(1,gpuDeviceCount);
        gpuMemory = zeros(1,gpuDeviceCount);
        gpuClock = zeros(1,gpuDeviceCount);
        gpuNumber = 1:gpuDeviceCount;
        for ii = 1:gpuDeviceCount
            gpu = gpuDevice(ii);
            gpuSupported(ii) = gpu.DeviceSupported;
            gpuMemory(ii) = gpu.TotalMemory / (1024^3);
            gpuClock(ii) = gpu.ClockRateKHz;
        end
        % sort out unsupported GPUs
        gpuNumber = gpuNumber(logical(gpuSupported));
        gpuMemory = gpuMemory(logical(gpuSupported));
        gpuClock = gpuMemory(logical(gpuSupported));
        % valuing memory over clock speed cause it usually is the limiting
        % factor
        if numel(find(gpuMemory == max(gpuMemory))) == 1
            % only one GPU has the highest amount of memory, use that one
            options.gpuSelection = gpuNumber(gpuMemory == max(gpuMemory));
        else
            % multiple GPUs have the highest amount of memory, find highest
            % clock speed
            gpuNumber = gpuNumber(gpuMemory == max(gpuMemory));
            gpuClock = gpuClock(gpuMemory == max(gpuMemory));
            options.gpuSelection = gpuNumber(gpuClock == max(gpuClock));
        end
        gpu = gpuDevice(options.gpuSelection);
        disp(['Selected ',gpu.Name])
        if gpu.TotalMemory/(1024^3) < 4
            warning('GPU has less than the recommended 4 GB of memory.')
        else
        end
    else
        % what even
        warning('Something''s odd with the GPU, falling back to CPU.')
        options.gpuEnabled = 'no';
    end
    
else
end

%% get pahts

filepath = paths.filepath;
savepath = paths.savepath;

% get base path of the file directory
% temporary results will be stored there, wrapper function will move the
% results to the appropriate location
pathParts = regexp(filepath,filesep,'split');
cutoff = numel(pathParts{end});
basePath = filepath(1:end-cutoff);



%% Single Evaluation

if exist(filepath,'file')
    % the file exists, check for segmentation
    [segments,segStatus] = ASAsegmentationWrapper(filepath,options);
    if segStatus ~= 0
        error('Error in file segmentation')
    else
    end
    if segments == 1
        evalstatus = ASAfullEvaluation(filepath,options);
        % move temp file to scan folder
        tempFile = fullfile(basePath,'EVALtemp.mat');
        destination = fullfile(savepath,'EVAL.mat');
        movefile(tempFile,destination)
        % generate report
        % not important for now, will be added later
        
        
        % plot all images and save them
        % l = load(destination);
        ASAfigures(destination,filepath,savepath,evalname,options)
        if strcmp(options.deleteEval,'yes')
            delete(destination)
        else
        end
    elseif segments > 1
        for jj = 1:segments
            segmentpath = [filepath(1:end-4),'-segment',num2str(jj),'.mat'];
            evalstatus = ASAfullEvaluation(segmentpath,options);
            % move temp file to scan folder
            savefolder = fullfile(savepath,['seg',num2str(jj)]);
            tempFile = fullfile(basePath,'EVALtemp.mat');
            mkdir(savefolder)
            destination = fullfile(savefolder,'EVAL.mat');
            movefile(tempFile,destination)
            % generate report
            % not important for now, will be added later
            
            
            % plot all images and save them
            % l = load(destination);
            ASAfigures(destination,segmentpath,savefolder,evalname,options)
            if strcmp(options.deleteEval,'yes')
                delete(destination)
            else
            end
        end
    else
        
    end
    % reset GPU if necessary
    if strcmp(options.gpuEnabled,'yes')
        reset(gpuDevice);
        gpu = gpuDevice(options.gpuSelection);
    else
    end
else
end



end



