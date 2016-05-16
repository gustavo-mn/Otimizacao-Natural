clear all
close all
clc

N_pop = 100; % Tamanho da população
L = 25; % Tamanho do genótipo
p_mut = 1/L; % Probabilidade de mutação
p_rec = 0.7; % Probabilidade de recombinação
N_ger = 100; % Número máximo de gerações
T = 10; % Quantidade de rodadas do algoritmo
abordagem = 'B'; % Abordagem de Baldwin ou de Lamarck na busca local
memetico = 1; % Decide se o algoritmo será memetico ou não

for t=1:T

    n = 1; % Geração atual
    fitness = zeros(N_ger,N_pop); % Vetor de fitness da população
    max_fitness = zeros(1,N_ger); % Vetor com a melhor fitness encontrada em cada geração
    min_fitness = zeros(1,N_ger); % Vetor com a pior fitness encontrada em cada geração
    media_fitness = zeros(1,N_ger); % Vetor com a média das fitness encontrada em cada geração

    P = unidrnd(2, [L, N_pop]) - 1; % Inicialização aleatória da população

    while (n <= N_ger)
        %% Cálculo de fitness
        
        fitness(n,:) = sum(P,1); % Cálculo das fitness de cada indivíduo
        
        max_fitness(n) = max(fitness(n,:));
        min_fitness(n) = min(fitness(n,:));
        media_fitness(n) = mean(fitness(n,:));

        if max_fitness(n) == 25 % Valor máximo de fitness é 25
            break;
        end
        
        if memetico % Aplica busca local
        
            fitness_melhorada = fitness(n,:); % Vetor de fitness após a busca local
            P_aux = P;
            
            switch (abordagem)
                case 'B' % Abordagem de Baldwin
                    for i = 1:size(P,1)
                        P_aux(i,:) = ~P_aux(i,:); % Inverte o valor dos bits correspondentes em cada indivíduo
                        
                        fitness_aux = sum(P_aux,1); % Calcula a fitness após modificação do bit i de cada indivíduo
                        
                        % A fitness de cada indivíduo é sempre a melhor encontrada até então, para cada um
                        
                        fitness_melhorada(fitness_aux > fitness_melhorada) = fitness_aux(fitness_aux > fitness_melhorada);
                        
                        P_aux = P; % Reseta a população para a original, de forma a alterar o próximo bit
                    end
                    
                case 'L' % Abordagem de Lamarck
                    P_melhores = P;
                    for i = 1:size(P,1)
                        P_aux(i,:) = ~P_aux(i,:); % Inverte o valor dos bits correspondentes em cada indivíduo
                        
                        fitness_aux = sum(P_aux,1); % Calcula a fitness após modificação do bit i de cada indivíduo
                        
                        % A fitness de cada indivíduo é sempre a melhor encontrada até então, para cada um
                        
                        fitness_melhorada(fitness_aux > fitness_melhorada) = fitness_aux(fitness_aux > fitness_melhorada);
                        
                        P_melhores(i, (fitness_aux > fitness_melhorada)) = P_aux(i, (fitness_aux > fitness_melhorada));
                        
                        P_aux = P; % Reseta a população para a original, de forma a alterar o próximo bit
                    end
                    
                    P = P_melhores; % Muda o genótipo para os indivíduos melhorados
            end
        end
        %% Seleção de pais

        % Probabilidade proporcional ao fitness
        
        if memetico
            pdf_fitness = fitness_melhorada/sum(fitness_melhorada);
            cdf_fitness = cumsum(pdf_fitness);
        else
            pdf_fitness = fitness(n,:)/sum(fitness(n,:));
            cdf_fitness = cumsum(pdf_fitness);
        end
        
        % Algoritmo SUS

        i = 1;
        membro_atual = i;
        r = unifrnd(0, 1/N_pop);    
        reprodutores = zeros(1,N_pop);

        while (membro_atual <= N_pop)
            while (r <= cdf_fitness(i))
                reprodutores(membro_atual) = i;
                r = r + 1/N_pop;
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
                P_new((c+1):end, (i+1)) = P((c+1):end, reprodutores(p1));

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

    if (n < N_ger)
        fitness = fitness(1:n,:);
        max_fitness = max_fitness(1:n);
        min_fitness = min_fitness(1:n);
        media_fitness = media_fitness(1:n);
        geracoes_otimas(t).n = n;
    else
        geracoes_otimas(t).n = N_ger;
    end
    
    geracoes_otimas(t).maximos = max_fitness;
    geracoes_otimas(t).minimos = min_fitness;
    geracoes_otimas(t).medias = media_fitness;
end

if memetico

    if exist(['dados_questao_2_' num2str(N_ger) '_geracoes_' abordagem '.mat'], 'file')
        delete(['dados_questao_2_' num2str(N_ger) '_geracoes_' abordagem '.mat'])
    end
    
    save(['dados_questao_2_' num2str(N_ger) '_geracoes_' abordagem '.mat'], 'geracoes_otimas');

else
    
    if exist('dados_questao_2_nao_memetico.mat','file')
        delete('dados_questao_2_nao_memetico.mat')
    end
    
    save('dados_questao_2_nao_memetico.mat', 'geracoes_otimas');
end
    