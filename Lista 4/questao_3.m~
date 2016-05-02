clear all
close all
clc

N_pop = 30; % Tamanho da população
n = 30; % Quantidade de dimensões da função de Ackley
N_exec = 1; % Quantidade de execuções independentes do algoritmo
N_filhos = 200; % Número de filhos criados por geração
N_avaliacoes = 30000; % Quantidade de vezes que a função de Ackley será avaliada
N_ger = ceil(N_avaliacoes / N_pop); % Quantidade de gerações
tau = 1/sqrt(n);
tau_linha = 1/sqrt(2*sqrt(n));
epsilon = 0.1; % Valor mínimo para o passo de mutação
min_valores = Inf(N_exec, 2);
%melhores_solucoes = zeros(2*n, N_avaliacoes, N_exec);

for t = 1:N_exec

    P = [unifrnd(-30,30,n,N_pop); epsilon*ones(n,N_pop)]; % Inicialização aleatória da população
    % contador = 0;
    ger = 1; % Geração atual
    opt = Inf;
    fitness = Inf(N_ger,N_pop);
    melhores_fitness = Inf(1, N_ger);
    media_fitness = zeros(1, N_ger);
    piores_fitness = Inf(1, N_ger);

    while (ger <= N_ger) && (min(fitness(ger,:)) ~= 0)
            
            % Cálculo dos fitness
            fitness(ger,:) = -20 * exp(-0.2 * sqrt((1/n) * sum(P(1:n,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P(1:n,:)), 1)) + 20 + exp(1);
            % contador = contador + N_pop;
            opt = find(fitness(ger,:) == min(fitness(ger,:)), 1);
            melhores_fitness(ger) = min(fitness(ger,:));
            media_fitness(ger) = mean(fitness(ger,:)); 
            piores_fitness(ger) = max(fitness(ger,:));

            if min_valores(t,1) > min(fitness(ger,:))
                min_valores(t,:) = [min(fitness(ger,:)) ger];
            end

            %min_valores(t, contador) = min(fitness);
            %melhores_solucoes(:,contador,t) = P(:,opt);

            % Mutação

            P_filhos = zeros(size(P,1), N_filhos); % 'N_pop' filhos de mutação e '(N_filhos - N_pop)' de recombinação

            %for i = 1:size(P, 2)
            
                % Mutação dos passos
                mutacoes = exp(tau * randn(n,N_pop)) .* (ones(n,1) * exp(tau_linha * randn(1,n))); % Alteração dos passos de mutação
                %mutacoes = P(31:end, i);
                %mutacoes = mutacoes .* (exp(tau_linha * randn(1)) * exp(tau * randn(n,1))); % Alteração dos passos de mutação
                %keyboard();
                if ~isempty(find((mutacoes < epsilon), 1)) % Se houver algum passo menor que o mínimo
                    mutacoes(mutacoes < epsilon) = epsilon; % Ajusta o passo para o mínimo aceitável
                end
                
                P_filhos((n+1):end, 1:N_pop) = mutacoes;

                % Mutação das soluções
                
                %P_filhos(:, i) = [P(1:30, i)+mutacoes; mutacoes];
                P_filhos(1:n, 1:N_pop) = P(1:n, 1:N_pop) + mutacoes .* randn(n,N_pop);

            %end

            % Recombinações

            for n_rec = 1:(N_filhos - N_pop)
                
                % Seleção dos pais                

                n_p1 = unidrnd(N_pop);
                n_p2 = unidrnd(N_pop);
                while (n_p2 == n_p1)
                    n_p2 = unidrnd(N_pop);
                end
                        
                sol_filho = zeros(n, 1);
                mut_filho = zeros(n, 1);
                sol_pai1 = P(1:n, n_p1);
                sol_pai2 = P(1:n, n_p2);
                mut_pai1 = P((n+1):end, n_p1);
                mut_pai2 = P((n+1):end, n_p2);

                % Recombinação discreta (soluções)
                
                z = unifrnd(0,1, n, 1);

                for i = 1:length(z)
                    if z(i) < 0.5
                        sol_filho(i) = sol_pai1(i);
                    else
                        sol_filho(i) = sol_pai2(i);
                    end
                end

                %sol_filho(z < 0.5) = sol_pai1(z < 0.5);
                %sol_filho(z >= 0.5) = sol_pai2(z >= 0.5);

                % Recombinação intermediária (passos)

                mut_filho = (mut_pai1 + mut_pai2)/2;

                P_filhos(:, (N_pop+n_rec)) = [sol_filho; mut_filho];
                    
            end

            % Seleção dos sobreviventes

            % Generacional com 'N_pop' filhos mais aptos

            fitness_filhos = ones(1,size(P_filhos,2));
            
            fitness_filhos = -20 * exp(-0.2 * sqrt((1/n) * sum(P_filhos(1:n,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P_filhos(1:n,:)), 1)) + 20 + exp(1);

            idx = zeros(1,N_pop);

            for i = 1:N_pop

                idx(i) = find(fitness_filhos == min(fitness_filhos), 1);
                fitness_filhos(idx(i)) = Inf;

            end

            P = P_filhos(:, idx);

            ger = ger + 1;

    end

    if (ger <= N_ger)
        fitness = fitness(1:ger, :);
    end

end

