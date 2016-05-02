N = 20;
% N = 1000000;

x = random('unif',-1,1, N, 1);

y = 4*sqrt(1-x.^2);

% y = 2 * (1./sqrt(1-x.^2)); % arccosseno

sum(y)/N % Média

% Resultado para N = 20
% ans = 3.2492

% Resultado para N = 1000000
% ans = 3.1415