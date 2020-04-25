function [segments,status] = ASAsegmentation(rawDataInput,metaData,filepath)

% try
    splitPath = regexp(filepath,filesep,'split');
    fileName = splitPath{end};
    fileBaseName = fileName(1:end-4);
    savepath = fullfile(splitPath{1:end-1});
    
    [ySize,~] = size(rawDataInput);
    shortSide = min(size(rawDataInput));
    if ySize == shortSide
        rawDataInput = rawDataInput';
        flipped = 1;
    else
        flipped = 0;
    end
    
    range = shortSide;
    
    [ySizeRest,~] = size(rawDataInput);
    segmentno = 0;
    
    rest = rawDataInput;
    
    while ySizeRest > 2*range
        scan = rest;
        segment = scan(1:range,1:end);
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
        segment = scan(end-(range-1):end,1:end);
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
        rest = scan(range+1:end-range-1,1:end);
        [ySizeRest,~] = size(rest);
    end
    
    if ySizeRest < range
        [ySize,~] = size(scan);
        midpoint = round(ySize/2);
        halfRange = round(range/2);
        segment = scan(midpoint-(halfRange-1):midpoint+(range-halfRange),1:end);
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
    elseif ySizeRest == range
        segment = rest;
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
    elseif ySizeRest > range
        segment = rest(1:range,1:end);
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
        segment = rest(end-(range-1):end,1:end);
        segmentno = segmentno+1;
        saveSubroutine(segment,segmentno,flipped,savepath,metaData,fileBaseName)
    end
    
    orderSegmentfiles(segmentno,savepath,fileBaseName)
    segments = segmentno;
    status = 0;
% catch ME
%     warning(['Encountered error in ASAsegmentation(): ',ME.msgtext])
%     status = -1;
%     return
% end

function orderSegmentfiles(segmentno,savepath,evaluationName)
    position = 'beginning';
    numbers = 1:segmentno;
    for ii = 1:segmentno
        oldName = fullfile(savepath,[evaluationName,'-segmentTEMP',num2str(ii),'.mat']);
        if strcmp(position,'beginning')
            segmentNumber = 0.5*ii + 0.5;
            position = 'end';
        elseif strcmp(position,'end')
            offset = -0.5*ii + 1;
            segmentNumber = numbers(end+offset);
            position = 'beginning';
        end
        newName = fullfile(savepath,[evaluationName,'-segment',num2str(segmentNumber),'.mat']);
        disp(oldName)
        disp(newName)
        movefile(oldName,newName)
    end
end

function saveSubroutine(segment,segmentno,flipped,savepath,metaData,evaluationName)
    if flipped
        segment = segment';
    else
    end
    rawData = segment;
    save(fullfile(savepath,[evaluationName,'-segmentTEMP',num2str(segmentno),'.mat']),...
        'rawData','metaData')
end
end