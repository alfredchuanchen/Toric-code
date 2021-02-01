load( 'parameters.mat' );
% ---import the parameterization of the model--- %

% index = 3;

thetas = theta_list( index, 1:4 );

h = H_BdG( N, t, Delta_epsilon, thetas(1), thetas(2), thetas(3), thetas(4) );

[ eigenvec, eigenval ] = eig( h );

[ eigenval, ind ] = sort( diag(eigenval), 'descend' );
eigenvec = eigenvec( :, ind );

% since we want the energies sligned in the form { E_1, E_2, ... , -E_1,
% -E_2, ... , -E_N }, with E_1 > E_2 > ... > E_N. We need to sort the
% negative energies in an ascending order:

[ negen, ind1 ] = sort( eigenval( N^2+1 : 2*N^2 ), 'ascend' );
eigenvec( :, N^2+1 : 2*N^2) = eigenvec( :, ind1 + N^2 );

eigenval( N^2+1 : 2*N^2 ) = negen;

% --done-- %


% ---define the new eigenvector matrix in the form of [ u, v^*; v, u^* ]:
% this is because the eigenvectors calculated numerically may differ by
% some phase factor like "-1"

neweigenvec = zeros( 2*N^2 );
neweigenvec( 1:2*N^2, 1:N^2 ) = eigenvec( 1:2*N^2, 1:N^2 );
neweigenvec( 1:N^2, N^2+1 : 2*N^2 ) = conj( eigenvec( N^2+1 : 2*N^2, 1:N^2 ) );
neweigenvec( N^2+1 : 2*N^2, N^2+1 : 2*N^2 ) = conj( eigenvec( 1 : N^2, 1:N^2 ) );

% ---export the data--- %

save( strcat( "spectra-N-",string(N),"-t-",string(t),...
    "-Delta_epsilon-",string(Delta_epsilon),"-index-",string(index),".dat" ),...
    "eigenval", '-ascii', '-double' );

save( strcat( "spectra-N-",string(N),"-t-",string(t),...
    "-Delta_epsilon-",string(Delta_epsilon),"-index-",string(index),".mat" ),...
    "neweigenvec" );

quit

% plot( zeros(size(eigenval)), eigenval, 'o' )