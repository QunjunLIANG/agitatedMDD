clear all
close all
clc

%% Setup
% Set parameters
num_nodes      = 214;    % number of time series, also, number of ROIs`
wkDir          = '/home/mri/MDD_kangning_projects/MDD_symptom_specific/inputs/timeseries_Power_typical/';
signal_file    = '_timeseries_power264.csv'; % input the postfix of the file name
outdir         = 'ETS_estimation';    % set directory for saving out images here
mkdir([wkDir outdir])

% Specify data parameters
sbjNamePath    = '/home/mri/MDD_kangning_projects/MDD_symptom_specific/inputs/';
subject_list   = importdata([sbjNamePath 'subject_list_for_TD.tsv']); %!!!!!!!!!!!!!!!!!!!!
subjects       = subject_list;

% how much subject you want to use for estimation, assign 1:numel(subj1ects) for alll subjects
useSubject     = 1:numel(subjects); 
outFile_prefix = 'ETS'; % add string to the output csv to identify the analysis
parallel_core  = 4; % number of core for parallel

%% working in loop 
p = parpool(parallel_core);
parfor s = useSubject
    tic
    subj = subjects{s};
    disp(['Processing ' subj]);
    
    % load regional time series
    ts = load([wkDir subj signal_file]);

	% z-score time series
	ts = zscore(ts);
	
    % get dimensions
    [ntime,nnodes] = size(ts);

    % calculate number of edges
    nedges = nnodes*(nnodes - 1)/2;

    % indices of unique edges (upper triangle)
    [u,v] = find(triu(ones(nnodes),1));
    idx = (v - 1)*nnodes + u;

    %% calculate static fc
    fc = corr(ts);

    %% generate edge time series
    ets = ts(:,u).*ts(:,v);

    %% calculate event amplitude

    % calculate co-fluctuation amplitude at each frame
    rms = sum(ets.^2,2).^0.5;

    % fraction of high- and low-amplitude frames to retain
    frackeep = 0.1;
    nkeep = round(ntime*frackeep);

    % sort co-fluctuation amplitude
    [~,idxsort] = sort(rms,'descend');

    % estimate fc using just high-amplitude frames
    fctop = corr(ts(idxsort(1:nkeep),:));

    % do the same using just low-amplitude
    fcbot = corr(ts(idxsort(end - nkeep + 1:end),:));

    fig(s) = figure('position',[560,528,700,420]);
    subplot(3,2,2); plot(rms);
    xlabel('frames'); ylabel('rss');

    subplot(3,2,1); imagesc(fc,[-1,1]); axis square;
    title('fc all time points');

    subplot(3,2,3); imagesc(fctop,[-1,1]); axis square;
    title('fc high-amplitude time points');

    subplot(3,2,4); scatter(fctop(idx),fc(idx),'k.'); xlim([-1,1]); ylim([-1,1]);
    text(-1,1,sprintf('r = %.2f',corr(fctop(idx),fc(idx))));

    subplot(3,2,5); imagesc(fcbot,[-1,1]); axis square;
    title('fc low-amplitude time points');

    subplot(3,2,6); scatter(fcbot(idx),fc(idx),'k.'); xlim([-1,1]); ylim([-1,1]);
    text(-1,1,sprintf('r = %.2f',corr(fcbot(idx),fc(idx))));
    
    %% export image
    saveas(fig(s), [wkDir outdir '/' subj '_ETS_plot.png'])
    
    %% export data
    csvwrite([wkDir outdir '/' subj '_rms.csv'], rms);
    csvwrite([wkDir outdir '/' subj '_FC.csv'], fc);
    csvwrite([wkDir outdir '/' subj '_bottomFC.csv'], fcbot);
    csvwrite([wkDir outdir '/' subj '_topFC.csv'], fctop);

end

