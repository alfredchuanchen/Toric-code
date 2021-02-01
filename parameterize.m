% ---parameter of the model--- %

N = 50;
t = 1;
Delta_epsilon = 10;

% ---generate the theta list--- %
L = 50; % ---steps for theta to change from 0 to pi.

nrows = 4*L-3;
theta_list = zeros( nrows, 4 );

theta_list( 1 : L, 1 ) = linspace( 0, pi, L );
theta_list( L+1 : nrows, 1 ) = pi;

theta_list( L : 2*L-1, 2 ) = linspace( 0, pi, L );
theta_list( 2*L : nrows, 2 ) = pi;

theta_list( 2*L-1 : 3*L-2, 3 ) = linspace( 0, pi, L );
theta_list( 3*L-1 : nrows, 3 ) = pi;

theta_list( 3*L-2 : nrows, 4 ) = linspace( 0, pi, L );

% ---done--- %

save( 'parameters.mat', 'theta_list', 'L', 'N', 't', 'Delta_epsilon' );