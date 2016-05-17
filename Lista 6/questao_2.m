clear all
close all
clc

N_pop = 100; % Tamanho da popula��o
L = 25; % Tamanho do gen�tipo
p_mut = 1/L; % Probabilidade de muta��o
p_rec = 0.7; % Probabilidade de recombina��o
N_ger = 100; % N�mero m�ximo de gera��es
T = 10; % Quantidade de rodadas do algoritmo
abordagem = 'B'; % Abordagem de Baldwin ou de Lamarck na busca local
memetico = 1; % Decide se o algoritmo ser� memetico ou n�o

for t=1:T

    n = 1; % Gera��o atual
    fitness = zeros(N_ger,N_pop); % Vetor de fitness da popula��o
    max_fitness = zeros(1,N_ger); % Vetor com a melhor fitness encontrada em cada gera��o
    min_fitness = zeros(1,N_ger); % Vetor com a pior fitness encontrada em cada gera��o
    media_fitness = zeros(1,N_ger); % Vetor com a m�dia das fitness encontrada em cada gera��o

    P = unidrnd(2, [L, N_pop]) - 1; % Inicializa��o aleat�ria da popula��o

    while (n <= N_ger)
        %% C�lculo de fitness
        
        fitness(n,:) = sum(P,1); % C�lculo das fitness de cada indiv�duo
        
        max_fitness(n) = max(fitness(n,:));
        min_fitness(n) = min(fitness(n,:));
        media_fitness(n) = mean(fitness(n,:));

        if max_fitness(n) == 25 % Valor m�ximo de fitness � 25
            break;
        end
        
        if memetico % Aplica busca local
        
            fitness_melhorada = fitness(n,:); % Vetor de fitness ap�s a busca local
            P_aux = P;
            
            switch (abordagem)
                case 'B' % Abordagem de Baldwin
                    for i = 1:size(P,1)
                        P_aux(i,:) = ~P_aux(i,:); % Inverte o valor dos bits correspondentes em cada indiv�duo
                        
                        fitness_aux = sum(P_aux,1); % Calcula a fitness ap�s modifica��o do bit i de cada indiv�duo
                        
                        % A fitness de cada indiv�duo � sempre a melhor encontrada at� ent�o, para cada um
                        
                        fitness_melhorada(fitness_aux > fitness_melhorada) = fitness_aux(fitness_aux > fitness_melhorada);
                        
                        P_aux = P; % Reseta a popula��o para a original, de forma a alterar o pr�ximo bit
                    end
                    
                case 'L' % Abordagem de Lamarck
                    P_melhores = P;
                    for i = 1:size(P,1)
                        P_aux(i,:) = ~P_aux(i,:); % Inverte o valor dos bits correspondentes em cada indiv�duo
                        
                        fitness_aux = sum(P_aux,1); % Calcula a fitness ap�s modifica��o do bit i de cada indiv�duo
                        
                        % A fitness de cada indiv�duo � sempre a melhor encontrada at� ent�o, para cada um
                        
                        fitness_melhorada(fitness_aux > fitness_melhorada) = fitness_aux(fitness_aux > fitness_melhorada);
                        
                        P_melhores(i, (fitness_aux > fitness_melhorada)) = P_aux(i, (fitness_aux > fitness_melhorada));
                        
                        P_aux = P; % Reseta a popula��o para a original, de forma a alterar o pr�ximo bit
                    end
                    
                    P = P_melhores; % Muda o gen�tipo para os indiv�duos melhorados
            end
        end
        %% Sele��o de pais

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

        % Reprodu��o

        P_new = zeros(size(P)); % Nova gera��o

        for i = 1:2:(size(P,2) - 1)
            
            p1 = unidrnd(length(reprodutores));
            p2 = unidrnd(length(reprodutores));

            while (p2 == p1)
                p2 = unidrnd(length(reprodutores)); % Evita que a mesma posi��o do vetor de reprodutores seja sorteada
            end

            r = unifrnd(0,1);
            
            if (r < p_rec) % Haver� recombina��o
                c = unidrnd(19); % Define o ponto de corte para recombina��o
                
                P_new(1:c, i) = P(1:c, reprodutores(p1));
                P_new((c+1):end, i) = P((c+1):end, reprodutores(p2));
                P_new(1:c, (i+1)) = P(1:c, reprodutores(p2));
                P_new((c+1):end, (i+1)) = P((c+1):end, reprodutores(p1));

            else % Os pais ser�o somente copiados para a gera��o seguinte
                
                P_new(:, i) = P(:, reprodutores(p1));
                P_new(:, (i+1)) = P(:, reprodutores(p2));

            end
        end

        % Muta��o bit a bit

        for j = 1:size(P_new, 2)
            for i = 1:size(P_new, 1)
                r = unifrnd(0,1);
                
                if (r < p_mut) % Haver� muta��o
                    if P_new(i,j) == 0
                        P_new(i,j) = 1;
                    else
                        P_new(i,j) = 0;
                    end
                end
            end
        end

        %% Sele��o dos sobreviventes    

        P = P_new; % Sele��o Generacional
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
    