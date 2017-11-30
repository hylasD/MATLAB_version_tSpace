function spdists = spdists_lknn( spdists, l, verbose )

	remove_edges = [];

	for idx = 1:length( spdists )
		% remove l-k neighbors at random
		neighs = find( spdists( :, idx ) );
		k = length( neighs ); % count number of neighbors
		remove_indices = neighs( randsample( length( neighs ), k - l ) );
		idx_remove_edges = sub2ind( size( spdists ), remove_indices, ones( k - l, 1 ) * idx );
		remove_edges = [ remove_edges; idx_remove_edges ];

		if( verbose )
			if( mod( idx, 50000 ) == 0 )
				fprintf( 1, '%3.2f%%', idx / length( spdists ) * 100 );
			elseif( mod( idx, 10000 ) == 0 )
				fprintf( 1, '.' );
			end
		end
	end

	spdists( remove_edges ) = 0;
end
