X = [0 4 6 9];

Y = [3 3.4]';

T = [1 0.1 50];

%% Valor inicial da função J

% Inicialização aleatória da matriz p(y|x)

pyx = rand(1,length(X));
aux = 1 - pyx;
pyx = [pyx; aux];

disp('Matriz p(y|x) inicial:')
pyx

% Cálculo da distorção D

D = 0;

for i = 1:length(X)
    D = D + sum(pyx(:,i) .* (X(i) - Y).^2);
end

D = D / length(X);

disp('Distorção D inicial:')
D

% Cálculo da entropia H

Hx = - log(1/length(X));

Hyx = 0;

for i = 1:length(X)
    Hyx = Hyx + sum(pyx(:,i) .* log(pyx(:,i)));
end

Hyx = - Hyx / length(X);

H = Hx + Hyx;

disp('Entropia H inicial:')
H

% Cálculo da função custo J

T_ini = T(1);

J = D - T_ini*H;

disp(['Valor inicial de J em T = ' num2str(T_ini) ':'])
J

%% Loop das temperaturas

for k = 1:length(T)

    %% Cálculo da matriz p(y|x) (cond. part.)
    
    pyx = zeros(length(Y), length(X));

    for c = 1:length(X)

        pyx (:,c) = exp(-((X(c) - Y).^2)/T(k));
        Zx = sum(pyx(:,c));
        pyx (:,c) = pyx(:,c)/Zx;

    end

    pyx

    %% Cálculo da distorção D
    
    D = 0;

    for i = 1:length(X)
        D = D + sum(pyx(:,i) .* (X(i) - Y).^2);
    end

    D = D / length(X);

    D

    
    %% Cálculo da entropia H
    
    Hx = - log(1/length(X));

    Hyx = 0;

    for i = 1:length(X)
        Hyx = Hyx + sum(pyx(:,i) .* log(pyx(:,i)));
    end

    Hyx = - Hyx / length(X);

    H = Hx + Hyx

    %% Cálculo da função custo J

    J = D - T(k)*H
    
    %% Cálculo dos centróides y (cond. centróide)
    
    Y_new = zeros(size(Y));

    for i = 1:length(Y)
        Y_new(i) = (pyx(i,:) * X') / sum(pyx(i,:));
    end

    Y_new
    
    %% Plot dos resultados
    
    h1 = plot(X, zeros(size(X)), 'bo ', 'LineWidth', 20); 
    hold on; 
    h2 = plot(Y_new, zeros(size(Y_new)), 'r+ ', 'LineWidth', 20); 
    plot(-2:10, zeros(1,13), 'k')
%     plot(X, zeros(size(X)), 'k'); 
    hold off
    grid on
    set(gca, 'ytick', [], 'FontSize', 30)
    legend([h1 h2], {'Amostras X', 'Centróides Y'}, 'FontSize', 30)
    title(['Disposição dos centróides Y em T = ' num2str(T(k))], 'FontSize', 30)
    
    pause

end

