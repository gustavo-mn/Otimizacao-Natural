X = zeros(1,4);
X(1) = 1; % Estado Inicial

M = [0.50 0.25 0.25; 0.25 0.50 0.25; 0.25 0.25 0.50];

p0 = [0 1 0]';
p = p0;

for i = 1:3
    p = M * p;
    r = rand();
    
    p_cdf = cumsum(p);
    
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
            break
        end
    end
    
    X(i+1) = x;
    
%     % Código para plotar a PDF e a CDF de p
    
    n = 0:2; 
    y = p'; 
    figure(1); 
    stem(n,y, 'LineWidth', 2); 
    axis([-1 3 0 1])
    set(gca, 'FontSize', 30, 'xtick', [0 1 2], 'ytick', [0.25 0.50 0.75 1]); 
    xlabel('x', 'FontSize', 30); 
    ylabel(['p' num2str(i)], 'FontSize',30); 
    title(['PDF de p em t = ' num2str(i)], 'FontSize', 30);
    figure(2); 
    stairs(n, cumsum(y), 'LineWidth', 2);
    set(gca, 'FontSize', 30, 'xtick', [0 1 2], 'ytick', [0.25 0.50 0.75 1]);
    xlabel('x', 'FontSize', 30);
    ylabel(['Soma cumulativa das probabilidades de p' num2str(i)], 'FontSize', 30)
    title(['CDF de p em t = ' num2str(i)], 'FontSize', 30);
    
    pause
    
    p
    r
    x
end

X