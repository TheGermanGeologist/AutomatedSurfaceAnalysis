function [segments,status] = ASAsegmentationWrapper(filepath,options)

% load file in question
l = load(filepath);
rawData = l.rawData;
metaData = l.metaData;

% check if segmentation is needed
[ySize,xSize] = size(rawData);
if ySize == xSize
    % already isometric, exiting
    segments = 1;
    status = 0;
else
    % rectangular, segment or crop
    switch options.autoSegmentation
        case 'yes'
            % Segmentation
            [segments,substatus] = ASAsegmentation(rawData,metaData,filepath);
            if substatus == 0
                disp(['Automatically segmented data in file ',filepath,' into ',num2str(segments),' isometric segments'])
                status = 0;
            else
                warning('Error in ASAsegmentation() in ASAsegmentationWrapper()')
                status = -1;
                return
            end
        case 'no'
            % Cropping
            if strcmp(options.autoCropping,'no')
                warning(['Data in file ',filepath,' is not isometric, but automatic segmentation and cropping are disabled. Falling back to auto-cropping.'])
            else
            end
            % Cropped
            shortSide = min(size(rawData));
            rawData = rawData(1:shortSide,1:shortSide);
            % backup old data, save new data
            backupName = [filepath(1:end-4),'original.mat'];
            movefile(filepath,backupName)
            save(filepath,'rawData','metaData')
            segments = 1;
            status = 0;
        otherwise
            warning('Invalid preference in: options.autoSegmentation. Must be ''yes'' or ''no''')
            status = -1;
            return
    end
end