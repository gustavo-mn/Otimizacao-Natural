X = [0 4 6 9];

n = 1;
D = zeros(1,12);

for t = -1:10
    
    i = find(X <= t);
    
    if isempty(i) || (length(X) == length(i))
        y1 = mean(X);
        y2 = y1;
    else
        y1 = sum(X(i))/length(i);
        y2 = sum(X((length(i)+1):end))/(length(X) - length(i));
    end
    
    D(n) = 0;
    contador = 1;
    
    for j = 1:length(X)
        if isempty(i)
            D(n) = D(n) + (X(j) - y1) ^ 2;
        else
            if contador > length(i)
                D(n) = D(n) + (X(j) - y2) ^ 2;
            else
                D(n) = D(n) + (X(j) - y1) ^ 2;
            end
            
            contador = contador + 1;
        end  
    end
    
    D(n) = (1/length(X)) * D(n);
    
    n = n + 1;
end

stem(-1:10, D);
set(gca, 'FontSize', 30, 'xtick', -1:10)
title('Valores de D(t)', 'FontSize', 30);
xlabel('t', 'FontSize', 30);
ylabel('D(t)', 'FontSize', 30);