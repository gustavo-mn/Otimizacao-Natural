Npop = 30; % Tamanho da população
n = 30; % Quantidade de dimensões da função de Ackley
Nexec = 1; % Quantidade de execuções independentes do algoritmo
Nfilhos = 200; % Número de filhos criados por geração
Navaliacoes = 200000; % Quantidade de vezes que a função de Ackley será avaliada
tau = 1/sqrt(n);
tau_linha = 1/sqrt(2*sqrt(n));
epsilon = 0.1; % Valor mínimo para o passo de mutação
min_valores = Inf(Nexec, 2);
%melhores_solucoes = zeros(2*n, Navaliacoes, Nexec);

for t = 1:Nexec

    P = unifrnd(-30,30,2*n,Npop); % Inicialização aleatória da população
    contador = 0;
    opt = Inf;
    fitness = Inf(1,Npop);

    while (contador < Navaliacoes) && (min(fitness) ~= 0)
            
            % Cálculo dos fitness
            fitness = -20 * exp(-0.2 * sqrt((1/n) * sum(P(1:30,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P(1:30,:)), 1)) + 20 + exp(1);
            contador = contador + 1;
            opt = find(fitness == min(fitness), 1); 

            if min_valores(t,1) > min(fitness)
                min_valores(t,:) = [min(fitness) contador];
            end

            %min_valores(t, contador) = min(fitness);
            %melhores_solucoes(:,contador,t) = P(:,opt);

            % Mutação

            P_filhos = zeros(size(P,1), Nfilhos); % 'Npop' filhos de mutação e '(Nfilhos - Npop)' de recombinação

            %for i = 1:size(P, 2)
            
                % Mutação dos passos
                mutacoes = exp(tau * randn(n,Npop)) .* (ones(n,1) * exp(tau_linha * randn(1,n))); % Alteração dos passos de mutação
                %mutacoes = P(31:end, i);
                %mutacoes = mutacoes .* (exp(tau_linha * randn(1)) * exp(tau * randn(n,1))); % Alteração dos passos de mutação
                %keyboard();
                if ~isempty(find((mutacoes < epsilon), 1)) % Se houver algum passo menor que o mínimo
                    mutacoes(mutacoes < epsilon) = epsilon; % Ajusta o passo para o mínimo aceitável
                end
                
                P_filhos(31:end, 1:Npop) = mutacoes;

                % Mutação das soluções
                
                %P_filhos(:, i) = [P(1:30, i)+mutacoes; mutacoes];
                P_filhos(1:30, 1:Npop) = P(1:30, 1:Npop) + mutacoes;

            %end

            % Recombinações

            for n_rec = 1:(Nfilhos - Npop)
                
                % Seleção dos pais                

                n_p1 = unidrnd(Npop);
                n_p2 = unidrnd(Npop);
                while (n_p2 == n_p1)
                    n_p2 = unidrnd(Npop);
                end
                        
                sol_filho = zeros(n, 1);
                mut_filho = zeros(n, 1);
                sol_pai1 = P(1:30, n_p1);
                sol_pai2 = P(1:30, n_p2);
                mut_pai1 = P(31:end, n_p1);
                mut_pai2 = P(31:end, n_p2);

                % Recombinação discreta (soluções)
                
                z = unifrnd(0,1, n, 1);

                sol_filho(z < 0.5) = sol_pai1(z < 0.5);
                sol_filho(z >= 0.5) = sol_pai2(z >= 0.5);

                % Recombinação intermediária (passos)

                mut_filho = (mut_pai1 + mut_pai2)/2;

                P_filhos(:, (Npop+n_rec)) = [sol_filho; mut_filho];
                    
            end

            % Seleção dos sobreviventes

            % Generacional com 'Npop' filhos mais aptos

            fitness_filhos = ones(1,size(P_filhos,2));
            
            fitness_filhos = -20 * exp(-0.2 * sqrt((1/n) * sum(P_filhos(1:30,:).^2, 1))) - exp((1/n) * sum(cos(2 * pi * P_filhos(1:30,:)), 1)) + 20 + exp(1);

            idx = zeros(1,Npop);

            for i = 1:Npop

                idx(i) = find(fitness_filhos == min(fitness_filhos), 1);
                fitness_filhos(idx(i)) = Inf;

            end

            P = P_filhos(:, idx);

    end
end

