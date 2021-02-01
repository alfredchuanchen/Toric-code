load('parameters.mat');

% index = 4;

index1 = index;

load( strcat( 'spectra-N-',string(N),'-t-',string(t),'-Delta_epsilon-',...
string(Delta_epsilon),'-index-',string(index1),'.mat' ) );

% ---notice here the directory of the eigenvector data: in the cluster, it
% is stored in the same directory as this script

mat1 = neweigenvec;
u1 = mat1( 1:N^2, 1:N^2 );
v1 = mat1( N^2+1 : 2*N^2, 1:N^2 );
clear neweigenvec;

if index1 < 4*L-3
    index2 = index1 + 1;

    load( strcat( 'spectra-N-',string(N),'-t-',string(t),'-Delta_epsilon-',...
        string(Delta_epsilon),'-index-',string(index2),'.mat' ) );

    mat2 = neweigenvec;
    u2 = mat2( 1:N^2, 1:N^2 );
    v2 = mat2( N^2+1 : 2*N^2, 1:N^2 );
    clear neweigenvec;

    M = [ u2'*conj(v2), -u2'*u1;...
        u1.'*conj(u2), v1.'*u1 ];
    M = (M-M.')/2; % -- here I am manually anti-symmetrizate M to same
    % time for the numerical check of whether M is skew-symmetric or not in
    % the pfaffian_LTL function
    
    op = (-1)^(N^2*(N^2+1)/2)*sqrt( abs(det(u1))*abs(det(u2)) )/...
        ( det(u1)*conj(det(u2)) )*pfaffian_LTL( M );
    % one should be careful about the power on top of '-1' in the formula
    % we used, should be the linear size of the matrix u, which is N^2,
    % i.e., system size!
else
    index2 = mod(index1 + 1, 4*L-3 );
    
    load( strcat( 'spectra-N-',string(N),'-t-',string(t),'-Delta_epsilon-',...
        string(Delta_epsilon),'-index-',string(index2),'.mat' ) );
    
    mat2 = neweigenvec;
    u2 = mat2( 1:N^2, 1:N^2 );
    v2 = mat2( N^2+1 : 2*N^2, 1:N^2 );
    clear neweigenvec;
    
    I_p = eye(N^2);
    I_p( funindex(N/2+1,N/2,N), funindex(N/2+1,N/2,N) ) = -1;
    
    newu2 = I_p*u2*I_p;
    newv2 = I_p*v2*I_p;
    clear  u2 v2;
    
    M = [ newu2'*conj(newv2), -newu2'*u1;...
        u1.'*conj(newu2), v1.'*u1 ];
    M = (M-M.')/2; % -- here I am again manually anti-symmetrizate M for the
    % same reason as before
    op = (-1)^(N^2*(N^2+1)/2)*sqrt( abs(det(u1))*abs(det(newu2)) )/...
        ( det(u1)*conj(det(newu2)) )*pfaffian_LTL( M );
    % here we use the LTL method to calculate the pfaffian
end

save( strcat( "overlap-N-",string(N),"-t-",string(t),...
    "-Delta_epsilon-",string(Delta_epsilon),"-index-",string(index1),".mat" ),...
    "op" );

quit


