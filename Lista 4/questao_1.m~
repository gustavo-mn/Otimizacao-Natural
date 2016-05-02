clear all
close all
clc

N_pop = 10; % Tamanho da população
N_bits = 20; % Quantidade de bits por genótipo
N_exec = 5; % Quantidade de vezes que o algoritmo será executado
p_rec = 0.8; % Probabilidade de haver recombinação
p_mut = 0.05; % Probabilidade de haver mutação
N_ger = 20; % Quantidade de gerações que serão consideradas
cont_ger = zeros(1,N_exec); % Armazena em que geração o ponto ótimo foi alcançado em cada execução

for t = 1:N_exec

    P = unidrnd(2, [N_bits, N_pop]) - 1; % Inicializa aleatóriamente a população
    n = 1; % Geração atual
    y = zeros(1, N_pop); % Valores da função y para cada indivíduo
    fitness = zeros(1, N_pop); % Fitness de cada indivíduo
    max_fitness = zeros(1, N_ger); % Maior valor de fitness encontrado até então na geração n
    min_fitness = zeros(1, N_ger); % Menor valor de fitness encontrado até então na geração n
    media_fitness = zeros(1, N_ger); % Média dos fitness na geração n

    while (n <= N_ger) % Critério de parada

        %% Cálculo da função y(x) para cada indivíduo

        for i = 1:size(P,2)
            
            G = P(:,i); % Genótipo do indivíduo
            ind = 1;
            contador = zeros(1,N_bits/4); % Conta quantos 1's existem nos segmentos do genótipo

            for j = 1:(N_bits/4)

                seg = G((4*(j-1)+1):4*j);
                contador(ind) = sum(seg);
                ind = ind + 1;

            end

            if (sum(contador) == 0) % Caso do genótipo ter todos os bits 0
                seg_vencedor = N_bits/4; % Prioridade dada ao segmento dos bits menos significativos
            else
                seg_vencedor = find(contador==max(contador)); % Encontra o segmento com maior contagem de 1's
                if (size(seg_vencedor,2) > 1) % Há segmentos empatados
                    seg_vencedor = min(seg_vencedor); % Prioridade para o segmento dos bits mais significativos
                end
            end

            switch(seg_vencedor) % Decodifica o genótipo em um fenótipo x
                case 1
                    x = -2;
                case 2
                    x = -1;
                case 3
                    x = 0;
                case 4
                    x = 1;
                case 5
                    x = 2;
            end
            
            y(i) = x^2 - 0.3*cos(10*pi*x); % Calcula o valor de y para o indíviduo
            fitness(i) = - 1/y(i); % Calcula o fitness do indivíduo
        end

        min_fitness(n) = min(fitness);
        max_fitness(n) = max(fitness);
        media_fitness(n) = mean(fitness);

        %% Seleção de pais

        % Probabilidade de cada indivíduo proporcional ao fitness

        pdf_fitness = abs(fitness)/sum(abs(fitness)); % Vetor de probabilidades de cada indivíduo

        cdf_fitness = cumsum(pdf_fitness);


        % Algoritmo 'Stochastic Universal Sampling'

        i = 1;
        membro_atual = i;
        r = unifrnd(0, 1/N_pop); % Sorteia um número da distribuição uniforme entre 0 e 1/N_pop
        reprodutores = zeros(1,N_pop); % Pais que serão selecionados para reproduzir

        while (membro_atual <= N_pop)
            while (r <= cdf_fitness(i))
                reprodutores(membro_atual) = i;
                r = r + 1/N_pop;
                membro_atual = membro_atual + 1;
            end
            i = i+1;
        end

        P_new = zeros(size(P)); % Nova geração

        for i = 1:2:(size(P_new)-1)


            p1 = unidrnd(length(reprodutores)); % Seleciona o primeiro pai
            p2 = unidrnd(length(reprodutores)); % Seleciona o segundo pai
            
            while (p2 == p1)
                p2 = unidrnd(length(reprodutores)); % Evita que seja sorteado o mesmo pai
            end
            
            % Recombinação uniforme
            
            r = unifrnd(0,1); % Sorteia um número da distribuição uniforme entre 0 e 1
            if (r < p_rec) % Ocorrerá recombinação
                for j = 1:N_bits
                    q = unifrnd(0,1); % Sorteia um número da distribuição uniforme entre 0 e 1
                    if (r > 0.5) 
                        P_new(j,i) = P(j, reprodutores(p1)); % Filho 1 recebe gene j do pai 1
                        P_new(j,(i+1)) = P(j, reprodutores(p2)); % Filho 2 recebe gene j do pai 2
                    else
                        P_new(j,i) = P(j, reprodutores(p2)); % Filho 1 recebe gene j do pai 2
                        P_new(j,(i+1)) = P(j, reprodutores(p1)); % Filho 2 recebe gene j do pai 1
                    end
                end
            else
                P_new(:,i) = P(:,reprodutores(p1));
                P_new(:,(i+1)) = P(:,reprodutores(p2));
            end
 
            % Recombinação de 1 crossover

           % r = unifrnd(0,1); % Sorteia um número da distribuição uniforme entre 0 e 1
           % if (r < p_rec) % Ocorrerá recombinação
           %     q = unidrnd(19);  % Sorteia um número entre 1 e 19, correspondente ao ponto de corte da recombinação
           %     P_new(1:q,i) = P(1:q,reprodutores(p1));
           %     P_new((q+1):end,i) = P((q+1):end,reprodutores(p2));
           %     P_new(1:q,(i+1)) = P(1:q,reprodutores(p2));
           %     P_new((q+1):end,(i+1)) = P((q+1):end,reprodutores(p1));
           % else
           %     P_new(:,i) = P(:,reprodutores(p1));
           %     P_new(:,(i+1)) = P(:,reprodutores(p2));
           % end

        end
        
        % Mutação bit a bit
        
        for j = 1:size(P_new,2)
            for k = 1:size(P_new,1)
                r = unifrnd(0,1); 
                if (r < p_mut) % Ocorrerá mutação
                    if (P_new(k,j) == 1)
                        P_new(k,j) = 0; % Comuta o bit
                    else
                        P_new(k,j) = 1; % Comuta o bit
                    end
                end
            end
        end

        % Seleção dos sobreviventes

        P = P_new; % Seleção Generacional
        n = n + 1;
    end

   % if (n < N_ger)
   %     max_fitness = max_fitness(1:(n-1));
   %     min_fitness = min_fitness(1:(n-1));
   %     media_fitness = media_fitness(1:(n-1));
   % end

    % Plotagem de gráficos
    
   % stem(max_fitness);
   % xlabel('Gerações', 'FontSize', 30);
   % ylabel('Valor máximo de fitness', 'FontSize', 30);
   % title('Valores máximos de fitness em cada geração', 'FontSize', 30);
   % set(gca, 'FontSize', 30);
   % figure;

   % stem(min_fitness);
   % xlabel('Gerações', 'FontSize', 30);
   % ylabel('Valor mínimo de fitness', 'FontSize', 30);
   % title('Valores mínimos de fitness em cada geração', 'FontSize', 30);
   % set(gca, 'FontSize', 30);
   % figure;

   % stem(media_fitness);
   % xlabel('Gerações', 'FontSize', 30);
   % ylabel('Valor médio de fitness', 'FontSize', 30);
   % title('Valores médios de fitness em cada geração', 'FontSize', 30);
   % set(gca, 'FontSize', 30);

   % keyboard();

   % close all

   cont_ger(t) = find(max_fitness == max(max_fitness),1);

end

