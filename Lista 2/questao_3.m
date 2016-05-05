T = zeros(1,20);
T0 = 0.1;
for i = 1:20
    T(i) = T0/log2(1+i);
end
% T = [0.100 0.0631 0.0500 0.0431 0.0387 0.0356 0.0333 0.0315 0.0301 0.0289];
x_atual = random('unif', -5,5,10,1);
% J_aux = sin(x_atual)./x_atual;
% J_aux(isnan(J_aux)) = 1;
% J_atual = sum(J_aux);

J_atual = sum(x_atual.^2);

J_min = J_atual;
N = 10000;
epsilon = 0.1;

J = zeros(length(T), N);
X = zeros(size(x_atual,1), N, length(T));

for k = 1:length(T)
    for n = 1:N
        
        r = random('unif',-1,1,size(x_atual));
        x_futuro = x_atual + epsilon * r;
        
%         J_aux = sin(x_futuro)./x_futuro;
%         J_aux(isnan(J_aux)) = 1;
%         J_futuro = sum(J_aux);
        
        J_futuro = sum(x_futuro.^2);
        
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
        
       J(k, n) = J_atual;
       X(:,n,k) = x_atual;
    end
    k       
end

% Resultados:
% 
% J_min =
% 
%     0.0075
%     
% X_min =
% 
%     0.0245
%     0.0190
%     0.0210
%    -0.0458
%     0.0146
%     0.0163
%    -0.0092
%    -0.0528
%    -0.0259
%     0.0022