clear all;
close all;
clc;

N_pop = 100; % Tamanho da população
N_ger = 20; % Quantidade de gerações
N_col = 8; % Número de colunas do tabuleiro
N_exec = 100; % Número de execuções independentes do algoritmo
p_rec = 0.8;
p_mut = 0.1;
ger_otimas = zeros(1,N_exec); % Armazena a geração em que a solução ótima ocorreu


for t = 1:N_exec

    fitness = zeros(N_ger,N_pop); % Número de rainhas salvas
    n = 1; % Geração atual
    max_fitness = 0;
    %% Inicialização aleatória dos genótipos

    P = zeros(N_pop, N_col);

    P(1,:) = 1:N_col;

    for i = 2:N_pop

        P(i,:) = P((i-1),:);

        r1 = unidrnd(N_col);
        r2 = unidrnd(N_col);
        while (r2 == r1)
            r2 = unidrnd(N_col);
        end
        
        aux = P(i,r1);
        P(i,r1) = P(i, r2);
        P(i,r2) = aux;
    end

    while (n <= N_ger) && (max_fitness ~= 8)

        %% Decodificação em fenótipo e cálculo de fitness

        F = zeros(N_col, N_col, N_pop);
        for i = 1:size(P,1);

            k = 1;
            for j = 1:size(P,2);
                
                F(P(i,j), k, i) = 1;
                k = k + 1;
                
            end
            
            k = 1;

            N_rainhas_salvas = 0;
            for j = 1:size(F,2)

                mask = zeros(N_col, N_col);
                mask(P(i,j), k) = 1;
                
                [lin,col] = find(mask == 1);         
                
                for caso = 1:4
                    lin_aux = lin;
                    col_aux = col;

                    switch(caso)
                        
                        case 1 % Diagonal superior direita
                        
                            while (lin_aux > 1) && (col_aux < N_col)
                                lin_aux = lin_aux - 1;
                                col_aux = col_aux + 1;
                                mask(lin_aux, col_aux) = 1;
                            end

                        case 2 % Diagonal superior esquerda
                
                            while (lin_aux > 1) && (col_aux > 1)
                                lin_aux = lin_aux - 1;
                                col_aux = col_aux - 1;
                                mask(lin_aux, col_aux) = 1;
                            end

                        case 3 % Diagonal inferior esquerda

                            while (lin_aux < N_col) && (col_aux > 1)
                                lin_aux = lin_aux + 1;
                                col_aux = col_aux - 1;
                                mask(lin_aux, col_aux) = 1;
                            end
                        
                        case 4 % Diagonal inferior direita

                            while (lin_aux < N_col) && (col_aux < N_col)
                                lin_aux = lin_aux + 1;
                                col_aux = col_aux + 1;
                                mask(lin_aux, col_aux) = 1;
                            end
                    end
                end
                if ((sum(sum(mask .* F(:,:,i))) - 1) == 0) % Sem rainhas na diagonal
                    N_rainhas_salvas = N_rainhas_salvas + 1;
                end
                
                k = k + 1;
            end

            fitness(n, i) = N_rainhas_salvas;

        end

        if max_fitness < max(fitness(n,:))
            max_fitness = max(fitness(n,:));
            opt = find(fitness(n,:) == max_fitness, 1);
            F_otimo = F(:,:,opt);
        end

        opt5 = find(fitness(n,:) == 5, 1);
        if ~isempty(opt5)
            F_3 = F(:,:,opt5);
        end

        %% Seleção de pais

        pdf_fitness = fitness(n,:) / sum(fitness(n,:));
        cdf_fitness = cumsum(pdf_fitness);

        % Algoritmo 'SUS'

        i = 1;
        membro_atual = i;
        r = unifrnd(0, 1/N_pop);
        reprodutores = zeros(1, N_pop);

        while (membro_atual <= N_pop)
            while (r <= cdf_fitness(i))
                reprodutores(membro_atual) = i;
                r = r + 1/N_pop;
                membro_atual = membro_atual + 1;
            end
            i = i+1;
        end

        % Recombinação

        P_filhos = zeros(size(P));
        
        for i = 1:2:(N_pop-1)

            n_p1 = unidrnd(N_pop);
            n_p2 = unidrnd(N_pop);
            while (n_p2 == n_p1) % Sorteia novamente o segundo pai
                n_p2 = unidrnd(N_pop);
            end

            p1 = P(reprodutores(n_p1), :);
            p2 = P(reprodutores(n_p2), :);

            r = unifrnd(0,1);
            if (r < p_rec) % Haverá recombinação

                % Crossover de 1 ponto            

                q = unidrnd(size(P_filhos, 2) - 2) + 1; % Sorteia um número entre 2 e 18
                P_filhos(i, 1:q) = p1(1:q);
                P_filhos((i+1), 1:q) = p2(1:q);

                contador_p1 = q; % Quantidade de genes copiados para o filho i
                contador_p2 = contador_p1; % Quantidade de genes copiados para o filho (i+1)
                idx = q + 1; % Próximo gene a ser analisado em cada pai

                while (contador_p1 ~= N_col) || (contador_p2 ~= N_col) % Ainda não preencheu completamente os filhos
                    if (contador_p1 ~= length(p1)) % Filho i não preenchido totalmente
                        if (isempty(find(P_filhos(i,:) == p2(idx), 1))) % Gene de permutação não repetido no filho i
                            contador_p1 = contador_p1 + 1;
                            P_filhos(i,contador_p1) = p2(idx); % Copia o gene do pai p2 para o filho i
                        end
                    end
                    if (contador_p2 ~= length(p2)) % Filho (i+1) não preenchido totalmente
                        if (isempty(find(P_filhos((i+1),:) == p1(idx), 1))) % Gene de permutação não repetido no filho (i+1)
                            contador_p2 = contador_p2 + 1;
                            P_filhos((i+1),contador_p2) = p1(idx); % Copia o gene do pai p1 para o filho (i+1)
                        end
                    end

                    idx = idx + 1;
                    
                    if (idx > N_col) % Percorreu todo genótipo de cada pai
                        idx = 1; % Volta para o começo do genótipo de cada pai
                    end
                end
                
            else % Copia os pais para os filhos; sem recombinação
                P_filhos(i,:) = p1;
                P_filhos((i+1),:) = p2;
            end
        end
        
        % Mutação
        
        r = unifrnd(0,1);
        for i = 1:size(P_filhos, 1) % Percorre cada indivíduo
            if (r < p_mut) % Haverá mutação no indivíduo
                p1 = unidrnd(N_col);
                p2 = unidrnd(N_col);
                while (p2 == p1) % Sorteia novamente p2
                    p2 = unidrnd(N_col);
                end
        
                % Troca as posições de dois elementos

                aux = P_filhos(i, p1); 
                P_filhos(i, p1) = P_filhos(i, p2);
                P_filhos(i, p2) = aux;
            end
        end

        %% Seleção dos sobreviventes

        P = P_filhos; % Generacional
                 
        n = n + 1;
    end

    ger_otimas(t) = n-1;

end

