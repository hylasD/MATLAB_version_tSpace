% tspace interface for calculation
% ddermadi@stanford.edu

clear
% get current path
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep];

% get all subdirectories
path_with_subdirectories = genpath(curr_path );
addpath( path_with_subdirectories );
savepath;
currentfolder = '';
allfiles = uipickfiles('num',[1 inf],'out','cell', 'FilterSpec', currentfolder, 'REFilter', '\.csv$');
if isequal(allfiles,0) ~=0
   return 
end
    
[~, ~, exts] = cellfun(@fileparts, allfiles, 'UniformOutput', false);
[uexts, IA, IC] = unique(exts);

for i=1:numel(uexts) %load csv file
    files = allfiles(IC==i);
    [path, ~, ext] = fileparts(files{1});
    %Read in header
    fids = cellfun(@(f) fopen(f, 'r'), files, 'UniformOutput', false);  %opening files
    csvheaders = cellfun(@(f) fgetl(f), fids, 'UniformOutput', false);   %readinf first line in file
    cellfun(@(f) fclose(f), fids, 'UniformOutput', false);  %closing files

    % Convert  header to cell array
    csvheaders = cellfun(@(h) regexp(h, '([^,]*)', 'tokens'), csvheaders, 'UniformOutput', false);
    csvheaders = cellfun(@(h) cat(2, h{:}), csvheaders, 'UniformOutput', false);
    
    % skip empty columns, index, time, beads, and cell or event length columns
    hdrs = csvheaders{1};
    omit = {'""', '', '"index"','index','"Time"','Time','"Cell_length"', 'Cell_length', '"Event_length"', 'Event_length', 'beads', '"beads"', 'Beads', '"Beads"'};
    %Indices of columns to be removed 
    Colindex = [find(ismember(hdrs,omit)), find(~cellfun(@isempty,regexpi(hdrs,'barcod')))];
    csvdats = cellfun(@(fname) csvread(fname, 1, 0), files, 'UniformOutput', false);
    csvdats = csvdats{i};
    
    % For multiple samples, to store that information
    if (input('If your file contains multiple samples labeled in specific column, press 1; otherwise press 0:'))
        prompt = sprintf('Type in the name of the column that labels samples: ');
        str = input(prompt, 's');
        SampleInd = find(~cellfun(@isempty,regexpi(hdrs,str)));
        if(find(SampleInd))
            disp('Column successfully found')
            Colindex = [Colindex, SampleInd];
            %Save Sample column
            Sample = csvdats(:,SampleInd);
        else
            disp('Check your column names and repeat loading of files')
        end   
    end
    
    
    csvdats(:,Colindex)=[];
    
    hdrs(:,Colindex)=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %transform data first, so far arcsin and logicle
    dataType = input('Press \n[1] for CyTOF/FACS Data or \n[2] scRNAseq: ');
    if (dataType ==1)
        if (input('To transform data, press 1 \n if you have transformed data, press 0: '))
        transType = input('Press \n[1] for CyTOF Data (arcsinh) or \n[2] FACS Data (logicle, be patient): ');
            if (transType ==1)
            %CyTOF data transformation with arcsin
                csvdats = asinh(csvdats./5);
            elseif (transType == 2)
            %FACS data transformation with LOGICLE transformation
            %It takes into account minimum value of each channel
                parfor (p = 1:size(csvdats,2),feature('numCores'))
                    w = (4.5-log10(262144/abs(min(csvdats(:,p)))))/2;
                    if w < 0
                        obj = [LogicleTransform(262144,0.5,4.5,0)];
                        csvdats(:,p) = obj.transform(csvdats(:,p));
                    else 
                        obj = [LogicleTransform(262144,w,4.5,0)];
                        csvdats(:,p) = obj.transform(csvdats(:,p));
                    end
                end 
            end
        end
        csvdats_rem = csvdats;
        hdrs_rem=hdrs;    
        removed = 0;
        for j = 1:numel(hdrs)
            column_name = hdrs{j};
            prompt = sprintf('Load %s [Y/n]: ', column_name);
            str = input(prompt, 's');
            if ~isempty(str) && ~strcmp(str,'y') && ~strcmp(str,'Y')
            % mark header as removed
            hdrs{j} = strcat(column_name, '_x');
            % remove
            hdrs_rem(:,j-removed) = [];
            csvdats_rem(:,j-removed) = [];
            removed = removed + 1;
            end
        end
    elseif (dataType == 2)
        if (input('To transform data, press 1 \n if you have transformed data, press 0: '))
                csvdats = log10(csvdats+1);
                scaleType = input('Press \n[1] for centering and scaling \n[2] only scaling usinf root-mean-square: ');
                if scaleType == 1
                    %line below is equivavelnt to R function: scale(data, center = T, scale = T))
                    csvdats = (csvdats - mean(csvdats,1))./std(csvdats,1);
                elseif scaleType == 2
                %Scale only by dividing values with root-mean-square equivavelnt to R function: scale(data, center = F, scale = T))
                    csvdats = csvdats./rms(csvdats,1);
                end
        end
    end
        csvdats_rem = csvdats;
        hdrs_rem=hdrs;    

    sessionData = zeros(0, size(csvdats_rem, 2));
    gates = cell(numel(csvdats),4);

    [~, csvname, ~] = fileparts(files{i}); 

    %-- add data to giant matrix
    currInd = size(sessionData, 1);
    sessionData(currInd+1:currInd+size(csvdats_rem,1), 1:size(csvdats_rem,2)) = csvdats_rem(:, :);
        
    gates{i, 1} = char(csvname);
    gates{i, 2} = currInd+1:currInd+size(csvdats_rem,1);
    gates{i, 3} = hdrs_rem;
    gates{i, 4} = files{i}; % opt cell column to hold filename
end
disp 'Files loaded';

% for PCA on tspaces
numFirstW = size(hdrs_rem, 2) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trajectory Algorithm Defaults
i = input('Select distance metric for tspace: [1]correlation (default), [2]cosine, [3]euclidean : ');
if i == 3
    parameters.metric = 'euclidean';
elseif i == 2
    parameters.metric = 'cosine';
else
    parameters.metric = 'correlation';  % default
end

parameters.k = input('Enter k number of neighbors (kNN): ');
if isempty(parameters.k)
    parameters.k = 20;
end

graph = input('Enter number of graphs: ');
if isempty(graph)
    graph = 5;
end

parameters.l = input('Enter l number of neighbors (l-kNN): ');
if isempty(parameters.l)
    L = parameters.k * 0.75;
    L = round(L, 0);
    parameters.l = L;
end

parameters.num_landmarks = input('Enter number of landmarks : ');
if isempty(parameters.num_landmarks)
    parameters.num_landmarks = 10;
end
parameters.partial_order = [];

voting = input('choose a voting scheme from the following: [1]uniform, [2]exponential (preferred), [3]linear, [4]quadratic : ');
if voting == 4
    parameters.voting_scheme = 'quadratic';
elseif voting == 3
    parameters.voting_scheme = 'linear';
elseif voting == 1
    parameters.voting_scheme = 'uniform';
else
    parameters.voting_scheme = 'exponential'; % default
end

% Landmarks are subsampled at equidistanced bands away from the start point.  
% Uses the points shortest path distance over the graph
parameters.band_sample = 1;

% Uses median 2 times
parameters.flock_landmarks = 2;

%Decide if you want to plot PCA embedding of tSPACE
plot = input('To plot tSPACE in PCA press 1 \n otherwise press 0: ');

% if 1 it will print messages, if 0 not.
parameters.verbose = 1;

%parameters.partial_order = []; % default is empty; Denis comment: in general this parameter is not used anywhere
parameters.deblur = 0;

% parameters.snn = input('enter number of shared nearest neighbors (default: 1) : ');
parameters.snn = 1;

% Searches for connected components if graph is discontious
parameters.search_connected_components = 1;

parameters.knn = []; % prebuilt on first run lNN graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters.exclude_points = [];
parameters.gates = gates;
parameters.gate_index = 1;
parameters.gateContext = gates{1, 2};

hdrs_rem2 = hdrs_rem;

%choose ground truth or # of kMeans populations
time_start1 = tic;
t_mode = input('For ground truth tSPACE press 1 \n For aproximate tSPACE press 0: ');
if (t_mode == 1)
    numPop = size(sessionData,1);
    clusters_trajectories = (1:1:size(sessionData,1))';
else % do kmeans
    sprintf(strcat('\nNumber of trajectories (default 200): '));
    numPop = input('\nInput number of trajectories to calculate; \n(recommendation: 100-200)');
    if isempty(numPop)
        numPop = 200;
    end
    rng(1); % For reproducibility
    clusters_trajectories = kmeans(sessionData, numPop, 'MaxIter', 10000); % 'Options', options);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part that calls trajectory script
    %separate cells by kMeans/SOM Population
    indexPops = zeros(numPop, numPop);
    pop = zeros(numPop, 1);
    for i = 1:size(clusters_trajectories,1)
        pop(clusters_trajectories(i)) = pop(clusters_trajectories(i)) + 1;
        indexPops(clusters_trajectories(i),pop(clusters_trajectories(i))) = i;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build kNN graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % print parameters
if parameters.verbose
    parameters
    disp 'building kNN graph';
end

% Build kNN graph
    tic;
    knn = parfor_spdists_knngraph( sessionData, parameters.k,...
        'distance', parameters.metric,...
        'chunk_size', 1000,... % TODO: parameterize and add opt for ppl without PC toolbox
        'verbose', parameters.verbose );
    if parameters.verbose, fprintf('kNN computed: %gs\n', toc); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (parameters.deblur)
        [i, j, s] = find(knn);
        % flock each data point to its knn median
        for ith=1:numel(i)
            data(ith, :) = median(data(j(i==ith), :)); 
        end
        if parameters.verbose, fprintf('re-computing kNN after data median filter\n'); tic; end
    	
        knn = parfor_spdists_knngraph( data, parameters.k, 'distance', parameters.metric, 'chunk_size', 1000, 'SNN', true, 'verbose', true);
        
        if parameters.verbose, fprintf('kNN re-computed after data median filter: %gs\n', toc); end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shared Nearest Neighbor
    if (parameters.snn~=0)
        if (parameters.verbose), fprintf('updating using jaccard...\n'); tic; end

        [j, i, s] = find(knn);
        nData = size(sessionData,1);
        rem = cell(1, nData);

        tic;
        % for each point 1..n
        parfor ci=1:nData
            
            % grab neighbors            
            from = (ci-1)*parameters.k+1;
            to = ci*parameters.k;
            i_inds = from:to;
            i_neighs = j(i_inds);
            
            % for each neighbor
            for i_ind=i_inds
                i_neigh=j(i_ind);
                
                % grab the neighbor's neighbors
                from = (i_neigh-1)*parameters.k+1;
                to = i_neigh*parameters.k;
                j_neighs = j(from:to);
%                 j_neighs = j(i==i_neigh);
                
                % if the shared neighbors are not many enough
                if sum(ismember(i_neighs, j_neighs)) < parameters.snn
                    
                    % add them to remove list
                    rem{ci} = [rem{ci} i_ind];
                end
            end
        end

        rem = cell2mat(rem);

        % remove relevant indices
        i(rem) = [];
        j(rem) = [];
        s(rem) = [];
        knn = sparse(j, i, s);
    
        if parameters.verbose, fprintf('jaccard computed: %gs\n', toc); end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   knn = spdists_undirected( knn ); 
   parameters.knn = knn;
   parameters.search_connected_components = true; 
   % pre-allocate memeory
   tspacem = zeros(size(sessionData,1), numPop);
   graph_panel = cell(graph,1);
    
    for graph_iter = 1:graph
    
	    if (parameters.k~=parameters.l)
	        lknn = spdists_lknn( parameters.knn, parameters.l, parameters.verbose );
	    else
	        lknn = parameters.knn;
	    end
            
        parfor (i = 1:numPop, feature('numCores'))
            %start event
            s = indexPops(i,1);
            tspacem(:,i) = runpathFinder(sessionData, lknn, parameters, s);
        end
        
        graph_panel{graph_iter} = tspacem;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tspacem = cat(3, graph_panel{:});
    tspacem = mean(tspacem,3);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %run pca on tspaces
    rng(1); % For reproducibility
    [tCoeff, tScore, tLatent, tTsquared, tExplained, tMu] = pca(tspacem,'NumComponents', 10);
   time_end = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          SAVING FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % file generation:
    % hdrs_rem corresponds to tspaceData (sessionData headers, W_n)
    % hdrs_rem2 corresponds to sessionData (sessionData headers only)
    % hdrs_rem3 corresponds to tScore (PCn_explain)

    % generating informative file name strings for csv and figure titles
    if (t_mode == 1)
    Clusters = 'sc';
    else
    Clusters = num2str(numPop);
    end
    
   %creating header for kMeans
    hdrs_rem3 = 'Cluster';
    
    %creating headers for wPCA
    for i = 1:size(tScore,2)
        pcaname = sprintf('tPC%d',i);
        pcaname = strcat(pcaname,'_',num2str(round(tExplained(i,1),2)));
        hdr = cellstr(pcaname);
        hdrs_rem3 = horzcat(hdrs_rem3, hdr);
    end
    
    %Headrs for tSPACE matrix file
    traj_hdrs = [];
    for i = 1:numPop
     tname = sprintf('T_%d',i);
        hdr = cellstr(tname);
        traj_hdrs = horzcat(traj_hdrs, hdr);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %file with only tspacem matrix (trajectories matrix) or tSPACE
    fileName = matlab.lang.makeValidName(['tsp' datestr(clock, 'mmddyy_HHMM') '_tSPACE_'   parameters.metric 'K' num2str(parameters.k) 'L' num2str(parameters.l) 'G' num2str(graph) 'LM' num2str(parameters.num_landmarks) 'T' Clusters csvname]);
    fileName = strcat(fileName, '.csv');
   
    %manipulate headers
    fid = fopen(fileName, 'w');
    fprintf(fid, '%s,', traj_hdrs{1,1:end-1});
    fprintf(fid, '%s\n', traj_hdrs{1,end});
    fclose(fid);

    dlmwrite(fileName, tspacem, '-append');
    fprintf('your tSPACE Matrix is exported as csv file %s\n', fileName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %file with allData (Clusters, pca on tspace, original data)
    index = zeros(size(sessionData,1),1);
    for i = 1:size(sessionData,1)
        index(i,1) = i;
    end
   
    if exist('Sample', 'var')
        if exist('trajectoriesOriginal')
            SampleO=Sample;
            Sample=Sample(~any(isnan(trajectoriesOriginal),2),:);
        end
        allData = cat(2,clusters_trajectories, tScore, csvdats, index, Sample);
        hdrs_rem3 = horzcat(hdrs_rem3,hdrs);
        hdrs_rem3 = horzcat(hdrs_rem3, cellstr('Index'));
        hdrs_rem3 = horzcat(hdrs_rem3, cellstr('Sample'));
          
    else ~exist('Sample', 'var')
        allData = cat(2, clusters_trajectories, tScore, csvdats, index);
        hdrs_rem3 = horzcat(hdrs_rem3, hdrs);
        hdrs_rem3 = horzcat(hdrs_rem3, cellstr('Index'));  
    end

    fileNameAll = matlab.lang.makeValidName(['tsp' datestr(clock, 'mmddyy_HHMM') '_PCA_'   parameters.metric 'K' num2str(parameters.k) 'L' num2str(parameters.l) 'G' num2str(graph) 'LM' num2str(parameters.num_landmarks) 'T' Clusters csvname]);
    fileNameAll = strcat(fileNameAll, '.csv');
    fid = fopen(fileNameAll, 'w');
    fprintf(fid, '%s,', hdrs_rem3{1,1:end-1});
    fprintf(fid, '%s\n', hdrs_rem3{1,end});
    fclose(fid);

    dlmwrite(fileNameAll, allData,'-append');
    fprintf('PCA components of tSPACE are exported as csv file %s\n', fileNameAll);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (plot == 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plotting PCAs
    scatter3(tScore(:,1), tScore(:,2), tScore(:,3),1,'k')

    plotName = matlab.lang.makeValidName(['tsp' datestr(clock, 'mmdd_HHMM') 'tSPACE_PCA'   parameters.metric 'K' num2str(parameters.k) 'L' num2str(parameters.l) 'G' num2str(graph) 'LM' num2str(parameters.num_landmarks) 'T' Clusters csvname]);
    title (plotName);
    
    lx = xlabel(hdrs_rem3(1,2)); % x-axis label
    set(lx,'rotation',15);
    ly = ylabel(hdrs_rem3(1,3));  % y-axis label
    set(ly,'rotation',-25);
    zlabel(hdrs_rem3(1,4));  % z-axis label
    savefig(plotName);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
