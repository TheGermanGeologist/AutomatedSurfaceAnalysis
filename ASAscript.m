% AUTOMATED SURFACE ANALYSIS (ASA)
% v.0.1

% by Matthias Alexander Dörfler
% University of Freiburg, 2019
% (c) Matthias Dörfler, 2019

% This software is distributed with a Creative Commons license
% (Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0))
% The terms and conditions given in the license.txt file apply. What
% follows is only a short summary:
% - You are allowed to share and adapt the licensed software
% - You must correctly attribute the original creator
% - You are not allowed to use the software commercially
% - You must distribute your modified version of the software under the
% same license as the original.

% The third-party functions distributed with this software are subject to
% the licenses given by their respective copyright holders. These licenses
% are summarized in the license.txt file and given individually with the
% third party functions. Where the license of a third party function
% disagrees with the license applicable to the software as a whole, the
% conditions given in the license for the third party function apply.



% Liability / warranty:
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


%%
clear
clc
close all
format compact

%% Options

% evaluation type can be 'single' or 'batch'. batch will require you to supply
% a txt file with a list of files to process
options.evaluationType = 'batch';
% for single evaluations, the name of the evaluation has to be specified
evalname = 'Sample01 Scan42';


% ROUGHNESS CALCULATION TYPES
% available: "raw" (no changes to data),
%            "highpassed" (used for microscale roughness by taking a lowpass
%                       FFT surface as reference)
%            "tilt" (removes tilt from data via 1st order plane fit)
%            "custom" (uses reference mode specified below)
options.calculationTypes = ["raw","highpassed","tilt","custom"];

% TILT
% tiltremoval mode - 'substract' or 'rotate'
options.tiltMode = 'substract';

% HISTOGRAMS
% enable histograms of the topography
options.histogram = 'yes';
% mode of histogram bin generation: 'fixed' or 'dynamic'
options.histogramBinMode = 'fixed';
% edges of histogram with fixed bins [µm]
options.histogramEdges = -1000:2:1000;

% FOURIER
% Enable Fourier Highpassing?
options.fourierHighpassing = 'yes';
% dynamically adjust cutoff frequency (using autocorrelation length)
options.dynamicFFTthreshold = 'yes';
% cutoff frequency in cycles/mm for the FFT lowpass reference surface
options.fftThreshold = 0.5;
% extended Fourier analyis (APSD,autocorrelation)
options.autocorrelationEnabled = 'yes';
options.APSDenabled = 'yes';
options.highpassedAPSD = 'yes';
% size of APSD angle bins (10° is recommended)
options.apsdBinsize = 10;

% REFERENCE SURFACE
% reference surface for standard roughness calculation can be 'zero' or
% 'mean'
options.RawReferenceMode = 'mean';
options.TiltReferenceMode = 'mean';
% Enable custom roughness calculation?
options.CustomRoughness = 'yes';
% options: 'zero', 'mean', 'surf1', ... 'surf5'
% mean value or n-th order fitted surface. 'surf1' would be the same as
% tilt corrected roughness
options.CustomReferenceMode = 'surf5';

% FIX DATA
% data fixing enabled?
% this will be performed and stored anyways if the fourier calculation is
% enabled
options.interpolationEnabled = 'yes';
% outlier removal enabled?
% this will use the matlab algorithm isoulier() to detect statistical
% outliers (possible measurement errors) and remove them.
options.outlierDetection = 'no';
% automatically crops data to square format
options.autoCropping = 'no';
% automatcially segment rectangular scans into square subsections?
% if this is not enabled, autoCropping will be forced
options.autoSegmentation = 'yes';
% Fix physical dimension? White Light Interferometry specific.
options.fixWLI = 'yes';

% SEMIVARIOGRAM
% compute Semivariogram of the data? Part of the spatial analysis. Very
% time intensive!
options.semivariogram2DEnabled = 'no';
% range of the semivariogram in pixel. Must be even integer!
options.semivariogramRange = 250;

% GPU COMPUTING
% attempt computation on the GPU, which might be faster
options.gpuEnabled = 'yes';

% MINIMUM DISK SPACE
% the minimum free disk space below which the evaluation will be aborted
% (in gigabyte)
options.minimumDisk = 5;

% SAVING
% Matlab offers different options for saving .mat files, some of which are
% faster, but result in (slightly) bigger file sizes
% '-v6' - fastest, no compression
% '-v7' - slowest, best compression
% '-v7.3' - somewhere in between but closer to -v7. adheres to HDF5 (.mat)
% support for true HDF5 will be added in a future release
options.saveType = '-v6';
% delete after saving images?
options.deleteEval = 'yes';
% image formats:
% 'png' 'jpeg' 'tiff' 'bmp' 'epsc' 'svg'
options.autoSaveImages = 'yes';
options.imageFormats = ["png"];

% ???
options.sillyness = 'yes';


%% Get paths (list with files to process and save location)

if strcmp(options.evaluationType,'single')
    % get file
    [filename,filepath] = uigetfile('.mat','Select the file to import');
    if isempty(filepath) | any(filepath == 0)
        return
    else
        filepath = [filepath,filename];
    end
    paths.filepath = filepath;
elseif strcmp(options.evaluationType,'batch')
    % get folder
    evalpath = uigetdir('Select the folder containing the files to evaluate');
    if isempty(evalpath) | (evalpath == 0)
        return
    else
    end
    paths.evalpath = evalpath;
else
    error('Invalid evaluation type. Must be ''single'' or ''multi''')
end
    
savepath = uigetdir('Select the folder to save the output into');
if isempty(savepath) | (savepath == 0)
    return
else
end
paths.savepath = savepath;


%% Make sure all needed functions are available
% add subdirectories of cd to search path (should hopefully contain the
% toolbox!)
addpath(genpath(cd));
% test for two needed additional toolboxes
if strcmp(options.gpuEnabled,'yes')
    try
        test = gpuArray(rand(10));
        clearvars test
    catch ME
        warning('Unable to use GPU. Make sure you have the Parallel Computing toolbox installed as well as a supported GPU.')
        disp(['Error message: ',ME.message'])
        warning('Continuing on CPU.')
        options.gpuEnabled = 'no';
    end
    else
end

if strcmp(options.APSDenabled,'yes')
    try
        test = repmat([0 0 1 0 0], 5,1);
        filled = imfill(logical(test),[3 1],4);
        clearvars test filled
    catch ME
        warning('Unable to use fcn imfill(). Make sure you have the Image Processing toolbox installed.')
        disp(['Error message: ',ME.message])
        warning('Disabling APSD calculation.')
        options.APSDenabled = 'no';
        options.highpassedAPSD = 'no';
    end
else
end


%% Perform Evaluation

switch options.evaluationType
    case 'single'
        status = ASAsingleWrapper(paths,evalname,options);
    case 'batch'
        status = ASAbatchWrapper(paths,options);
    otherwise
        
end
