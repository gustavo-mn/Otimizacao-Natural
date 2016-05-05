N = 100;
T = 4;
count = zeros(T,3);
pins = zeros(T,3);


X = zeros(N,T);
M = [0.50 0.25 0.25; 0.25 0.50 0.25; 0.25 0.25 0.50];
p0 = [0.33 0.33 0.33]';

for n = 1:N
    
    p = p0;

    for i = 0:(T-1)
        
        if (i ~= 0) % Não é o estado inicial
            p = M * p;
            r = rand();

            p_cdf = cumsum(p); % CDF da distribuição p

            for j = 1:length(p_cdf)
                if r <= p_cdf(j)
                    switch j
                        case 1 
                            x = 0;
                        case 2
                            x = 1;
                        case 3
                            x = 2;
                    end
                break;
                end
            end
            
            X(n, i+1) = x;
            
        else % Estado inicial; apenas inicializa X(n, 1)
            X(n, i+1) = random('unid', 3) - 1; % Estado Inicial    
        end
    end
end

% Trecho necessário para a questão d

for t = 1:size(X,2)
    figure(t);
    [count(t,:), pins(t,:)] = hist(X(:,t), [0 1 2]);
    hist(X(:,t), [0 1 2])
    axis([-1 3 0 50])
    set(gca, 'FontSize', 30)
    title(['Histograma dos estados X para t = ' num2str(t-1)], 'FontSize', 30)
end
pause
for t = 1:size(X,2)
    figure(t);
    stem(pins(t,:),count(t,:)/sum(count(t,:)), 'LineWidth', 2)
    axis([-1 3 0 1])
    set(gca, 'FontSize', 30)
    title(['PDF dos estados X para t = ' num2str(t-1)], 'FontSize', 30)
end
