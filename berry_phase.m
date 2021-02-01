load('parameters.mat');

phase_list = zeros(4*L-3,1);

for i = 1 : 4*L-3
    load( strcat( '/Users/chuanchen/Dropbox/Chuan-Works/Toric code/',...
        'fermion-flux/chern-number/data/berry_phase/overlap-N-',string(N),...
        '-t-',string(t),'-Delta_epsilon-',string(Delta_epsilon),...
    '-index-',string(i),'.mat' ) );
    phase_list(i,1) = op;
    clear op;
end

z = prod(phase_list);

bp = -imag( log( z/abs(z) ) );



