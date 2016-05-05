T = [0.100 0.0631 0.0500 0.0431 0.0387 0.0356 0.0333 0.0315 0.0301 0.0289];
J = [0.5 0.2 0.3 0.1 0.4];
x_atual = random('unid',5); % Estado inicial aleatório
J_atual = J(x_atual);
J_min = Inf;
N = 1000;
X = zeros(length(T), N);
custos = zeros(length(T), N);

for k = 1:length(T)
    for n = 1:N
        x_futuro = random('unid',5);
        J_futuro = J(x_futuro);
        delta_J = J_futuro - J_atual;

        if delta_J < 0
            x_atual = x_futuro;
            J_atual = J_futuro;
        else
            r = rand();
            if r < exp(-(delta_J)/T(k))
                x_atual = x_futuro;
                J_atual = J_futuro;
            end
        end

        if J_atual < J_min
            J_min = J_atual;
            x_min = x_atual;
        end
        
        X(k, n) = x_atual;
        custos(k, n) = J_atual;
    end
end

for k = 1:length(T)
    figure(k)
    hist(X(k, (N-100):end), [0 1 2 3 4 5 6])
    axis([0 6 0 100])
    set(gca, 'FontSize', 30)
    title(['Histograma dos 100 últimos estados X para T = ' num2str(T(k))], 'FontSize', 30)
end