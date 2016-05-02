T = zeros(1,10);
T0 = 1;
for i = 1:10
%     T(i) = T0/log2(1+i);
    T(i) = T0/i;
end

N = 1000000;
epsilon = 0.05;
x_atual = random('unif', -2,2,20,1);

f1 = sum((x_atual.^2));
f2 = sum((x_atual.*(x_atual-2)+1.05));
f3 = sum((x_atual.*(x_atual+3)+(1.55)^2));
f4 = sum((x_atual.*(x_atual-4)+4.05));
f5 = sum((x_atual.*(x_atual+4)+4.05));

J_atual = f1 .* f2 .* f3 .* f4 .* f5;
J_min = J_atual;

J = zeros(length(T), N);
X = zeros(size(x_atual,1),N,length(T));

for k = 1:length(T)
    for n = 1:N
        
        r = trnd(1,size(x_atual));
%         r = random('unif',-1,1,size(x_atual));
        
        x_futuro = x_atual + epsilon * r;
        
        f1 = sum((x_futuro.^2));
        f2 = sum((x_futuro.*(x_futuro-2)+1.05));
        f3 = sum((x_futuro.*(x_futuro+3)+(1.55)^2));
        f4 = sum((x_futuro.*(x_futuro-4)+4.05));
        f5 = sum((x_futuro.*(x_futuro+4)+4.05));
        
        J_futuro = f1 .* f2 .* f3 .* f4 .* f5;
        
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

            