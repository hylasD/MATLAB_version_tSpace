%ddermadi@stanford.edu

function [G] = pathFinder(data, lknn, parameters, s)
%    G = pathFinder( data, parameters) 
%
%    parameters: structure for specifying algorithm parameters which can contain
%    zero or more of the following fields:
%      
%      [k]                      : size of neighborhood k closest points
%      [s]                      : index to the starting point in the data
%      
%      [verbose]                : messages are printed to Matlab's sto
%      [metric]                 : string - distance metric for constructing 
%                                   the nearest neighbor graphs
%      [snn]                    : shared nearest neighbor
%      [search_connected_components] : search for connected components if
%                                   graph is discontious
%      [knn]                    : precomputed knn
%      
%      
% 

% set up return structure
G.T = []; % traj
G.traj = [];
G.path = [];


rng('shuffle');
    
	% run traj. landmarks
	[ traj, dist, iter_l, paths_l2l] = trajectory_waypoint_p(data, lknn, G, parameters, s);
    
    G.landmarks = iter_l; % iter_l is l from trajectory_landmarks_mod function
    G.traj = {traj};
    G.dist = {dist};
    G.path = {paths_l2l};
	
    % calculate weighed trajectory
    if strcmpi(parameters.voting_scheme, 'uniform')
        T_full(:, :) = ones(numel(iter_l), size(data,1));
    elseif strcmpi(parameters.voting_scheme, 'exponential')
        sdv = mean ( std ( dist) )*3;
        T_full = exp( -.5 * (dist / sdv).^2);
    elseif strcmpi(parameters.voting_scheme, 'linear')
        T_full = repmat(max(dist), size( dist, 1 ), 1) - dist;
        if ~isempty(parameters.exclude_points)
            T_full(:, parameters.exclude_points) = 1;
        end
    elseif strcmpi(parameters.voting_scheme, 'quadratic')
        T_full = repmat(max(dist.^2), size( dist, 1 ), 1) - dist.^2;
    end
    
    % The weghing matrix must be a column stochastic operator
    T_full = T_full ./ repmat( sum( T_full ), size( T_full, 1 ), 1 );
        
    T = T_full;
    
    % save initial solution - start point's shortest path distances
    t=[];
    t( 1,:)  = traj(1,:);
	t( end+1, : ) = sum( traj .* T );
    
	% iteratively realign trajectory (because landmarks moved)
	converged = 0; user_break = 0; realign_iter = 2;

	while  ~converged && ~user_break
		realign_iter = realign_iter + 1;

		traj = dist;
        for idx = 1:size( dist, 1 )
			% find position of landmark in previous iteration
			idx_val = t( realign_iter - 1, iter_l( idx ) );
			% convert all cells before starting point to the negative
			before_indices = find( t( realign_iter - 1, : ) < idx_val );
			traj( idx, before_indices ) = -dist( idx, before_indices );
			% set zero to position of starting point
			traj( idx, : ) = traj( idx, : ) + idx_val;
        end
               
		% calculate weighed trajectory
		t( realign_iter, : ) = sum( traj .* T );

		% check for convergence
        fpoint_corr = corr( t( realign_iter, : )', t( realign_iter - 1, : )' );
        fprintf( 1, '%2.5f...', fpoint_corr);
		converged = fpoint_corr > 0.9999;
        
        if (mod(realign_iter,16)==0)
            % break after too many alignments - something is wrong
            user_break = true;
            %fprintf('\nWarning: Force exit after %g iterations\n', realign_iter);
        end
	end
    
%     plot_iterations(parameters.plot_data, t);
    
	fprintf( 1, '\n%d realignment iterations, ', realign_iter-1 );

	% save final trajectory for this graph    
    G.T = t(realign_iter, :);
    
end

