T = zeros(1,10);
T0 = 0.1;
for i = 1:10
    T(i) = T0/log2(1+i);
end

N = 10000;
epsilon = 0.1;

x_atual = random('unif', -5,5,10,1);

J_aux = sin(x_atual)./x_atual + (2*sin(x_atual - 10)./(x_atual - 10)).^10; % 'Limita' a sync(x) para valores de x entre -10 e 10
J_aux(isnan(J_aux)) = 1;
J_atual = sum(J_aux);

J_min = J_atual;

J = zeros(length(T), N);
X = zeros(size(x_atual,1),N,length(T));

for k = 1:length(T)
    for n = 1:N
        
        r = random('unif', -1,1, size(x_atual));
        
        x_futuro = x_atual + epsilon * r;
        J_aux = sin(x_futuro)./x_futuro + (2*sin(x_futuro - 10)./(x_futuro - 10)).^10; 
        J_aux(isnan(J_aux)) = 1;
        J_futuro = sum(J_aux);
        
        delta_J = J_futuro - J_atual;
        
        if delta_J < 0
            x_atual = x_futuro;
            J_atual = J_futuro;
        else
            a = rand();
            
            if a < exp(-(delta_J)/T(k))
                x_atual = x_futuro;
                J_atual = J_futuro;
            end
        end
        
        if J_atual < J_min
            J_min = J_atual;
            X_min = x_atual;
        end
        
        J(k,n) = J_atual;
        X(:,n,k) = x_atual;
        
    end
end

% Resultados
% 
% J_min =
% 
%    -2.1542
%    
% X_min =
% 
%    -4.3440
%     4.3660
%     4.2649
%    -4.5078
%    -4.5456
%     4.5105
%     4.6312
%    -4.3770
%     4.5351
%    -4.6876   

            