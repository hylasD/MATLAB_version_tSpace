function [ traj, dist, l, paths_l2l, diffdists ] = trajectory_waypoint_p(data, lknn,  G, parameters, s)

    nl = parameters.num_landmarks;

	if( length( nl ) == 1 )
        [dists, paths, ~] = graphshortestpath(lknn, s,'METHOD','Dijkstra', 'directed', true);
        
        % if not given landmarks list, decide on random landmarks
        n_opts = 1:size(data,1);
        if (parameters.band_sample)
            n_opts = [];
            window_size = .1;
            max_dist = max(dists);
            for prc = .998:-window_size:.08
                band = find(dists>=(prc-window_size)*max_dist & dists <=prc*max_dist); % & num_jumps_arr >= floor((prc-.05)*max_jumps) & num_jumps_arr <= ceil(prc*max_jumps));
                n_opts = [n_opts randsample( band, min(length(band), nl - 1 - length(parameters.partial_order)), false )];
            end
        end
        nl = randsample( n_opts, nl - 1 - length(parameters.partial_order) );
        
        % flock landmarks 
        if (parameters.flock_landmarks > 0)
        for k=1:parameters.flock_landmarks
            [IDX, ~] = knnsearch(data, data(nl, :), 'distance', parameters.metric, 'K', 20);     
            for i=1:numel(nl)
                nl(i) = knnsearch(data, median(data(IDX(i, :), :)), 'distance', parameters.metric); 
            end
        end
        end
    end

    diffdists = zeros(length(nl), length(nl));

    partial_order = [s;parameters.partial_order(:)]; % partial_order includes start point
	l = [ partial_order; nl(:) ]; % add extra landmarks if user specified
    
    % calculate all shortest paths
    paths_l2l = cell(length(l),1);
    for li = 1:length( l )
        [dist( li, : ), paths, ~] = graphshortestpath( lknn, l( li ),'METHOD','Dijkstra', 'directed', false );
        if sum(cellfun(@(x)isempty(x), paths(l))) 
            fprintf('\nWarning: found empty path');
        end
        paths_l2l(li) = {paths(l)};
        unreachable = (dist(li,:)==inf);
        unreachable(parameters.exclude_points) = 0;

        while (any(unreachable) && parameters.search_connected_components)
            fprintf(['\n Warning: %g were unreachable. try increasing l'...
                'or k.Your data is possibly non continous, ie '...
                'has a completely separate cluster of points.'...
                'Wanderlust will roughly estimate their distance for now \n'],...
                sum(unreachable));
            if (parameters.plot_debug_branch)
                figure('Color',[1 1 1]);
                scatter(parameters.plot_data(:,1), parameters.plot_data(:,2), 150,'.b');
                hold on;
                scatter(parameters.plot_data(l(li),1), parameters.plot_data(l(li),2), 150,'.g');
                scatter(parameters.plot_data(unreachable,1), parameters.plot_data(unreachable,2), 150,'.r');
                title('Unreachable in red (from gre');
            end
            % find closest unreachable point to reachable points.
            % connect it on the spdists. continue iteratively.
            unreachablei = find(unreachable);
            reachablei = find(~unreachable);
            cou = 0;
            while ~isempty(unreachablei)
                cou = cou+1;
                [idx, d] = knnsearch(data(unreachablei, :), data(reachablei, :));
                closest_reachable = d==min(d);
                
                %add connection to spdists
                lknn(reachablei(closest_reachable),...
                    unreachablei(idx(closest_reachable))) = min(d);
                lknn(unreachablei(idx(closest_reachable)),...
                    reachablei(closest_reachable)) = min(d);
                % move points from unreachable list to reachable
                reachablei(end+1:end+length(find(closest_reachable))) = ...
                    unreachablei(idx(closest_reachable));
                unreachablei(idx(closest_reachable)) = [];
                
                if ~mod(cou, 10)
                    break;
                end
            end
            [dist( li, : ), paths, ~] = graphshortestpath( lknn, l( li ),'METHOD','Dijkstra', 'directed', false );
            paths_l2l(li) = {paths(l)};
            unreachable = (dist(li,:)==inf);
        end
        
        if( parameters.verbose )
            fprintf( 1, '.' );
        end
    end
    
    if ~isempty(parameters.exclude_points)
        dist(:, parameters.exclude_points) = mean(mean(dist~=inf));
    end
    
    if any(any(dist==inf))
        dist(dist==inf) = max(max(dist~=inf));
        if (parameters.verbose)
            fprintf('\nwarning: some points remained unreachable (dist==inf)');
        end
    end
    
    % adjust paths according to partial order by redirecting
    nPartialOrder = length(partial_order);
    for radius = 1:nPartialOrder 
        for landmark_row = 1:nPartialOrder
            if (landmark_row + radius <= nPartialOrder)
                a = landmark_row;
                b = landmark_row + (radius-1);
                c = landmark_row + radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
            if (landmark_row - radius >= 1)
                a = landmark_row;
                b = landmark_row - (radius-1);
                c = landmark_row - radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
        end
    end

    
    
	% align to dist_1 - this for loop refers to partial order stuff %DOes
	% not do anything to regular wandelrust calcualtion
	traj = dist;
    for idx = 2:length(partial_order)
        [~, closest_landmark_row] = min(dist); %closest landmark will determine directionality
        traj(idx, closest_landmark_row < idx) = -dist(idx, closest_landmark_row < idx);
        traj( idx, : ) = traj( idx, : ) + dist( 1, l( idx ) );
    end
    
    % This is the actual align for regular wanderlust
    if length( l ) > length(partial_order)
        for idx = length(partial_order)+1:length( l )
            % find position of landmark in dist_1
            idx_val = dist( 1, l( idx ) );
            % convert all cells before starting point to the negative
            before_indices = find( dist( 1, : ) < idx_val );
            traj( idx, before_indices ) = -dist( idx, before_indices );
            % set zero to position of starting point
            traj( idx, : ) = traj( idx, : ) + idx_val;
        end
    end
%     if (parameters.plot_landmark_paths)
% %         plot_landmark_paths(parameters.plot_data, paths_l2l, l);
%     end
%     if (parameters.branch)
%         [RNK, bp, diffdists, Y] = splittobranches(traj, traj(1, :), data, l, dist, paths_l2l, parameters);
%     end
end