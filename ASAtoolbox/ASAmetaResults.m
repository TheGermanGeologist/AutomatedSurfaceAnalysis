function status = ASAmetaResults(evalfilepath,metafilepath,list,fileNo,segNo,options)

status = -1;

warning('off','MATLAB:table:RowsAddedExistingVars')

l = load(evalfilepath);

if exist(metafilepath,'file')
    % load table
    m = load(metafilepath);
    metaTable = m.metaResults;
    metaResults = fillASAtable(metaTable,l,list,fileNo,segNo,options);
    save(metafilepath,'metaResults');
else
    metaTable = initASAtable();
    metaResults = fillASAtable(metaTable,l,list,fileNo,segNo,options);
    save(metafilepath,'metaResults');
end

end

function filledTable = fillASAtable(metaTable,results,list,fileNo,segNo,options)
    pos = numel(metaTable.SampleName) + 1;
    % Sample / Scan information
    metaTable.SampleName{pos} = list.SampleName{fileNo};
    metaTable.ScanNo(pos) = list.ScanNo(fileNo);
    metaTable.SegNo(pos) = segNo;
    metaTable.Material{pos} = list.Material{fileNo};
    metaTable.MaterialPrefix{pos} = list.MaterialPrefix{fileNo};
    metaTable.TestType{pos} = list.TestType{fileNo};
    metaTable.StrainRate(pos) = list.Strainrate(fileNo);
    metaTable.StrainRateDev(pos) = list.StrainrateDev(fileNo);
    metaTable.sigma1(pos) = list.sigma1(fileNo);
    metaTable.sigma3(pos) = list.sigma3(fileNo);
    metaTable.size{pos} = size(results.topography);
    if isa(list.EMod,'cell')
        metaTable.EMod(pos) = str2double(list.EMod{fileNo});
    else
        metaTable.EMod(pos) = list.EMod(fileNo);
    end
    metaTable.ScanMethod{pos} = list.ScanMethod{fileNo};
    % Results
    if isfield(results,'percentOutl')
        metaTable.OutlierPercent(pos) = results.percentOutl;
        metaTable.NaNpercent(pos) = results.percentNaN;
    else
    end
    if any(strcmp(options.calculationTypes,'raw'))
        metaTable.RawRoughness{pos} = results.rawRoughness;
        if strcmp(options.histogram,'yes')
            metaTable.Nraw{pos} = results.Nraw;
            metaTable.edgesRaw{pos} = results.edgesRaw;
        else
        end
    else
    end
    if any(strcmp(options.calculationTypes,'tilt'))
        metaTable.TiltRoughness{pos} = results.tiltRoughness;
        if strcmp(options.histogram,'yes')
            metaTable.Ntilt{pos} = results.Ntilt;
            metaTable.edgesTilt{pos} = results.edgesTilt;
        else
        end
    else
    end
    if any(strcmp(options.calculationTypes,'highpassed'))
        metaTable.HighpassedRoughness{pos} = results.highpassedRoughness;
        if strcmp(options.histogram,'yes')
            metaTable.Nhighpassed{pos} = results.Nhighpassed;
            metaTable.edgesHighpassed{pos} = results.edgesHighpassed;
        else
        end
    else
    end
    if any(strcmp(options.calculationTypes,'custom'))
        metaTable.CustomRoughness{pos} = results.customRoughness;
        if strcmp(options.histogram,'yes')
            metaTable.Ncustom{pos} = results.Ncustom;
            metaTable.edgesCustom{pos} = results.edgesCustom;
            metaTable.CustomType{pos} = options.CustomReferenceMode;
        else
        end
    else
    end
    
    if strcmp(options.autocorrelationEnabled,'yes')
        metaTable.betaMin(pos) = results.betaMin(1);
        metaTable.betaMinAngle(pos) = results.betaMinAngle(1);
    else
    end
    
    if strcmp(options.APSDenabled,'yes')
        [~,dirPos] = max(results.APSD);
        metaTable.maxAPSDdir(pos) = results.binAngles(dirPos);
    else
    end
    
    if strcmp(options.highpassedAPSD,'yes')
        [~,dirPos] = max(results.highpassedAPSD);
        metaTable.maxHighpassedAPSDdir(pos) = results.binAngles(dirPos);
    else
    end
    
    if strcmp(options.semivariogram2DEnabled,'yes')
        metaTable.semiMeanDir(pos) = mean(results.truedir);
        metaTable.semiMeanEccent(pos) = mean(results.eccentr);
    else
    end
    
    filledTable = metaTable;
end

function newTable = initASAtable()
    % initialize table
    SampleName = {}; ScanNo = []; SegNo = []; Material = {}; MaterialPrefix = {};
    TestType = {}; StrainRate = []; StrainRateDev = []; sigma1 = []; sigma3 = [];
    EMod = []; ScanMethod = {}; RawRoughness = {}; TiltRoughness = {}; HighpassedRoughness = {};
    CustomRoughness = {}; CustomType = {}; Nraw = {}; edgesRaw = {}; Ntilt = {};
    edgesTilt = {}; Nhighpassed = {}; edgesHighpassed = {}; Ncustom = {};
    edgesCustom = {}; betaMin = []; betaMinAngle = []; maxAPSDdir = [];
    maxHighpassedAPSDdir = []; semiMeanDir = []; semiMeanEccent = [];
    size = {}; OutlierPercent = []; NaNpercent = [];
    
    newTable = table(SampleName,ScanNo,SegNo,Material,MaterialPrefix,TestType,...
        StrainRate,StrainRateDev,sigma1,sigma3,EMod,ScanMethod,size,...
        OutlierPercent,NaNpercent,RawRoughness,...
        TiltRoughness,HighpassedRoughness,CustomRoughness,CustomType,Nraw,...
        edgesRaw,Ntilt,edgesTilt,Nhighpassed,edgesHighpassed,Ncustom,edgesCustom,...
        betaMin,betaMinAngle,maxAPSDdir,maxHighpassedAPSDdir,semiMeanDir,semiMeanEccent);
end