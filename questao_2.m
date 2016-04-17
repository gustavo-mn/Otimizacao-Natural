Npop = 100; % Tamanho da população
L = 25; % Tamanho do genótipo
p_mut = 1/L; % Probabilidade de mutação
p_rec = 0.7; % Probabilidade de recombinação
Nger = 100; % Número máximo de gerações
n = 1; % Geração atual
T = 10; % Quantidade de rodadas do algoritmo

geracoes_otimas = zeros(1,T); % Guarda em que gerações alcançou-se o máximo de fitness

for t=1:T

    fitness = zeros(1,Npop); % Vetor de fitness da população
    max_fitness = zeros(1,Nger); % Vetor com a melhor fitness encontrada em cada geração
    min_fitness = zeros(1,Nger); % Vetor com a pior fitness encontrada em cada geração
    media_fitness = zeros(1,Nger); % Vetor com a média das fitness encontrada em cada geração

    P = unidrnd(2, [L, Npop]) - 1; % Inicialização aleatória da população

    while (n <= Nger) && (max_fitness(n) ~= 25) % Valor máximo de fitness é 25
        
        %% Cálculo de fitness
        
        fitness(n,:) = sum(P,1); % Cálculo das fitness de cada indivíduo

        max_fitness(n) = max(fitness(n,:));
        min_fitness(n) = min(fitness(n,:));
        media_fitness(n) = mean(fitness(n,:));

        %% Seleção de pais

        % Probabilidade proporcional ao fitness

        pdf_fitness = fitness/sum(fitness);
        cdf_fitness = cumsum(pdf_fitness);
        
        % Algoritmo SUS

        i = 1;
        membro_atual = i;
        r = unifrnd(0, 1/Npop);    
        reprodutores = zeros(1,Npop);

        while (membro_atual <= Npop)
            while (r <= cdf_fitness(i))
                reprodutores(membro_atual) = i;
                r = r + 1/Npop;
                membro_atual = membro_atual + 1; 
            end
            i = i + 1;
        end

        % Reprodução

        P_new = zeros(size(P)); % Nova geração

        for i = 1:2:(size(P,2) - 1)
            
            p1 = unidrnd(length(reprodutores));
            p2 = unidrnd(length(reprodutores));

            while (p2 == p1)
                p2 = unidrnd(length(reprodutores)); % Evita que a mesma posição do vetor de reprodutores seja sorteada
            end

            r = unifrnd(0,1);
            
            if (r < p_rec) % Haverá recombinação
                c = unidrnd(19); % Define o ponto de corte para recombinação
                
                P_new(1:c, i) = P(1:c, reprodutores(p1));
                P_new((c+1):end, i) = P((c+1):end, reprodutores(p2));
                P_new(1:c, (i+1)) = P(1:c, reprodutores(p2));
                P_new((c+1):end, (i+1)) = P((c+1):end, reprodutores(p2));

            else % Os pais serão somente copiados para a geração seguinte
                
                P_new(:, i) = P(:, reprodutores(p1));
                P_new(:, (i+1)) = P(:, reprodutores(p2));

            end
        end

        % Mutação bit a bit

        for j = 1:size(P_new, 2)
            for i = 1:size(P_new, 1)
                r = unifrnd(0,1);
                
                if (r < p_mut) % Haverá mutação
                    if P_new(i,j) == 0
                        P_new(i,j) = 1;
                    else
                        P_new(i,j) = 0;
                    end
                end
            end
        end

        %% Seleção dos sobreviventes    

        P = P_new; % Seleção Generacional
        n = n + 1;
    end

    if (n < Nger)
        max_fitness = max_fitness(1:n);
        min_fitness = min_fitness(1:n);
        media_fitness = media_fitness(1:n);
    end
    
    if (t == 1)
        % Plot dos gráficos

        stem(max_fitness);
        xlabel('Geração', 'FontSize', 30);
        ylabel('Fitness Máxima', 'FontSize', 30);
        title('Melhores fitness por geração', 'FontSize', 30);
        set(gca, 'FontSize', 30);
        figure;

        stem(min_fitness);
        xlabel('Geração', 'FontSize', 30);
        ylabel('Fitness Mínima', 'FontSize', 30);
        title('Piores fitness por geração', 'FontSize', 30);
        set(gca, 'FontSize', 30);
        figure;

        stem(media_fitness);
        xlabel('Geração', 'FontSize', 30);
        ylabel('Média Fitness', 'FontSize', 30);
        title('Média dos fitness por geração', 'FontSize', 30);
        set(gca, 'FontSize', 30);
    end

    geracoes_otimas(t) = n-1;

end

disp(['Média de gerações para alcançar o máximo: ' num2str(mean(geracoes_otimas))])
disp(['Desvio-padrão de gerações para alcançar o máximo: ' num2str(std(geracoes_otimas, 1))])
