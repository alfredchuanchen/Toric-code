function op = H_BdG(N,t,del_epsilon,theta1,theta2,theta3,theta4)


% ---define the matrix for the c^{\dagger}c hopping part--- %

Zeta = del_epsilon*eye(N^2);

% ---first define the hopping matrix for terms without the last row and column--- %

for y = 1 : N-1
    for x = 1 : N-1
        Zeta( funindex(x,y,N), funindex(x+1,y,N) ) = -t;
        Zeta( funindex(x+1,y,N), funindex(x,y,N) ) = -t;
        Zeta( funindex(x,y,N), funindex(x,y+1,N) ) = -t;
        Zeta( funindex(x,y+1,N), funindex(x,y,N) ) = -t;
    end
end

% ---including the last row and column--- %

for x = 1 : N-1
    Zeta( funindex(x,N,N), funindex(x+1,N,N) ) = -t;
    Zeta( funindex(x+1,N,N), funindex(x,N,N) ) = -t;
end

for y = 1 : N-1
    Zeta( funindex(N,y,N), funindex(N,y+1,N) ) = -t;
    Zeta( funindex(N,y+1,N), funindex(N,y,N) ) = -t;
end

% ---ed--- %

% ---branch cut to the left, flux is at the centre--- %

for x = 1 : N/2
    Zeta( funindex(x,N/2,N), funindex(x,N/2+1,N) ) = t;
    Zeta( funindex(x,N/2+1,N), funindex(x,N/2,N) ) = t;
end

Zeta( funindex(N/2+1,N/2,N), funindex(N/2+1,N/2+1,N) ) = -t*exp(1i*theta1);
Zeta( funindex(N/2+1,N/2+1,N), funindex(N/2+1,N/2,N) ) = -t*exp(-1i*theta1);

Zeta( funindex(N/2+1,N/2,N), funindex(N/2+2,N/2,N) ) = -t*exp(1i*theta2);
Zeta( funindex(N/2+2,N/2,N), funindex(N/2+1,N/2,N) ) = -t*exp(-1i*theta2);

Zeta( funindex(N/2+1,N/2,N), funindex(N/2+1,N/2-1,N) ) = -t*exp(1i*theta3);
Zeta( funindex(N/2+1,N/2-1,N), funindex(N/2+1,N/2,N) ) = -t*exp(-1i*theta3);

Zeta( funindex(N/2+1,N/2,N), funindex(N/2,N/2,N) ) = -t*exp(1i*theta4);
Zeta( funindex(N/2,N/2,N), funindex(N/2+1,N/2,N) ) = -t*exp(-1i*theta4);



% ---ed--- %



% ---define the pairing matrix for c^{dagger} c^{\dagger}--- %

Delta = zeros(N^2);

% ---first define the hopping matrix for terms without the last row and column--- %

for y = 1 : N-1
    for x = 1 : N-1
        Delta( funindex(x+1,y,N), funindex(x,y,N) ) = -t;
        Delta( funindex(x,y,N), funindex(x+1,y,N) ) = t;
        Delta( funindex(x,y,N), funindex(x,y+1,N) ) = -t;
        Delta( funindex(x,y+1,N), funindex(x,y,N) ) = t;
    end
end

% ---including the last row and column--- %

for x = 1 : N-1
    Delta( funindex(x+1,N,N), funindex(x,N,N) ) = -t;
    Delta( funindex(x,N,N), funindex(x+1,N,N) ) = t;
end

for y = 1 : N-1
    Delta( funindex(N,y,N), funindex(N,y+1,N) ) = -t;
    Delta( funindex(N,y+1,N), funindex(N,y,N) ) = t;
end

% ---ed--- %


% ---branch cut to the left, flux is at the centre--- %

for x = 1 : N/2
    Delta( funindex(x,N/2,N), funindex(x,N/2+1,N) ) = t;
    Delta( funindex(x,N/2+1,N), funindex(x,N/2,N) ) = -t;
end

Delta( funindex(N/2+1,N/2,N), funindex(N/2+1,N/2+1,N) ) = -t*exp(1i*theta1);
Delta( funindex(N/2+1,N/2+1,N), funindex(N/2+1,N/2,N) ) = t*exp(1i*theta1);

Delta( funindex(N/2+1,N/2,N), funindex(N/2+2,N/2,N) ) = t*exp(1i*theta2);
Delta( funindex(N/2+2,N/2,N), funindex(N/2+1,N/2,N) ) = -t*exp(1i*theta2);

Delta( funindex(N/2+1,N/2,N), funindex(N/2+1,N/2-1,N) ) = t*exp(1i*theta3);
Delta( funindex(N/2+1,N/2-1,N), funindex(N/2+1,N/2,N) ) = -t*exp(1i*theta3);

Delta( funindex(N/2,N/2,N), funindex(N/2+1,N/2,N) ) = t*exp(1i*theta4);
Delta( funindex(N/2+1,N/2,N), funindex(N/2,N/2,N) ) = -t*exp(1i*theta4);


% ---ed--- %


op = [ Zeta, Delta; Delta', -Zeta.' ];

end