clear all; close all; clc

N_pop = 200; % Tamanho da população
N_avaliacoes = 200000; % Quantidade de avaliações da função de Ackley
N_ger = N_avaliacoes/N_pop; % Número de gerações avaliadas 
n = 30; % Dimensões da função de Ackley
epsilon = 0.02; % Tamanho mínimo do passo de alteração
eta = 6; % Passo de alteração inicial
N_exec = 100; % Número de execuções do algoritmo
alpha = 0.2; % Constante de alteração do passo
q = 10; % Quantidade de vezes que uma solução participará do torneio
melhor_fitness_execucao = ones(1, N_exec); % Melhores fitness encontradas em cada execução do algoritmo
geracao_otima_execucao = zeros(1, N_exec); % Primeira geração em que apareceu a melhor fitness em cada execução do algoritmo

for t = 1:N_exec
    
    ger = 1; % Geração atual
    melhores_fitness = ones(1, N_ger); % Melhores fitness por geração 
    P = [unifrnd(-30,30,n,N_pop); eta*ones(n,N_pop)]; % Inicialização aleatória da população
    
    % Cálculo do fitness
    
    fitness = -20 * exp(-0.2 * sqrt((1/n) * sum(P(1:n,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P(1:n,:)), 1)) + 20 + exp(1); % Fitness por indivíduos 
    
    melhores_fitness(ger) = min(fitness);
    
    while (ger <= N_ger) && (melhores_fitness(ger) > 1e-7)
        
        % Mutação
        
        P_filhos = zeros(size(P));
        
        passos = P((n+1):end, :);
        passos = passos .* (1 + alpha * (ones(n, 1) * randn(1, N_pop))); % Mutação dos passos;
        passos(passos < epsilon) = epsilon;
        
        P_filhos((n+1):end, :) = passos;
        
        P_filhos(1:n, :) = P(1:n, :) + passos .* randn(n, N_pop);
        
        % Seleção de sobreviventes
        
        P_torneio = [P P_filhos]; % Monta a população de filhos + pai para o torneio
        fitness = -20 * exp(-0.2 * sqrt((1/n) * sum(P_torneio(1:n,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P_torneio(1:n,:)), 1)) + 20 + exp(1); % Fitness dos participantes
        resultados = zeros(1,size(P_torneio, 2)); % Resultados do torneio;
        qtd_participacoes = zeros(1,size(P_torneio, 2)); % Contagem de participações de cada indivíduo no torneio

        for p1 = 1:length(resultados)
            while (qtd_participacoes(p1) ~= q)
                
                qtd_participacoes(p1) = qtd_participacoes(p1) + 1;
                
                competidores_disponiveis = find(qtd_participacoes ~= q);                    
                p2 = unidrnd(size(P_torneio, 2));
                
                if length(competidores_disponiveis) > 1 % Há mais de um competidor disponível
                    while (p2 == p1) || (qtd_participacoes(p2) == q) % Evita sortear o mesmo participante ou um participante que já competiu o número máximo de vezes no torneio
                        p2 = unidrnd(size(P_torneio, 2));
                    end
                                        
                    qtd_participacoes(p2) = qtd_participacoes(p2) + 1;
                    
                    if fitness(p1) < fitness(p2)
                        resultados(p1) = resultados(p1) + 3;
                    else if fitness(p1) == fitness(p2)
                            resultados(p1) = resultados(p1) + 1;
                            resultados(p2) = resultados(p2) + 1;
                        else
                            resultados(p2) = resultados(p2) + 3;
                        end
                    end
                    
                else if length(competidores_disponiveis) == 1 % Sobrou somente 1 competidor; aceita sortear outro que já competiu 'q' vezes, porém não o pontua
                    while (p2 == p1) % Sorteia outro competidor, que não o próprio
                        p2 = unidrnd(size(P_torneio, 2));
                    end
                    
                    if fitness(p1) < fitness(p2)
                        resultados(p1) = resultados(p1) + 3;
                    else if fitness(p1) == fitness(p2)
                            resultados(p1) = resultados(p1) + 1;
                            resultados(p2) = resultados(p2) + 1;
                        end
                    end
                    end
                end                        
            end
        end
        
        idx = zeros(1, N_pop); % Índice que identifica os indivíduos melhor classificados
        
        for i = 1:N_pop % Salva os N_pop indivíduos melhores classificados no torneio para a próxima geração
            idx(i) = find(resultados == max(resultados), 1);
            P(:,i) = P_torneio(:,idx(i));
            resultados(idx(i)) = -Inf;
        end
        
        ger = ger + 1;
        
        fitness = fitness(idx);
        melhores_fitness(ger) = min(fitness);              
        
    end
    
    melhor_fitness_execucao(t) = min(melhores_fitness);
    geracao_otima_execucao(t) = find(melhores_fitness == min(melhores_fitness), 1);
end