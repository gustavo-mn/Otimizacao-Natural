x_atual = 0;
epsilon = 0.1;

T0 = 1;
N = 100000;
x = zeros(1,N);
J = zeros(1,N);
J_min = inf;
X_min = inf;

X = [];
J_rodadas = [];
k = 0;
interrompe = 0;
indice_T = 1;

while (~interrompe)
    
    T_atual = (0.9^k) * T0;
    
    disp(['T: ' num2str(T_atual)])
    
    if T_atual < 0.1
        interrompe = 1;
    end
    
    T(indice_T) = T_atual;
    indice_T = indice_T + 1;
    
    if (~interrompe)
    
        for i = 1:N

            x(i) = x_atual;
            J_atual = (((x_atual + 2) ^ 2) - 3) * (10 * ((x_atual - 3) ^ 2) + 5);
            J(i) = J_atual;

            if J_atual < J_min
                J_min = J_atual;
                X_min = x_atual;
            end

            r = random('unif', -1, 1);
            x_futuro = x_atual + epsilon * r;
            J_futuro = (((x_futuro + 2) ^ 2) - 3) * (10 * ((x_futuro - 3) ^ 2) + 5);
            J_dif = J_futuro - J_atual;

            if J_dif < 0
                x_atual = x_futuro;
            else
                beta = rand(1,1);

                if beta < exp(-(J_dif)/T_atual)
                    x_atual = x_futuro;
                end
            end        
        end

        X = [X; x];
        J_rodadas = [J_rodadas; J];

        k = k + 1;
    end
    
end

% Resultados do código
% 
% J_min = -845.6411
% X_min = -2.4936