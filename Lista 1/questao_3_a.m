x_atual = random('norm', 0, sqrt(0.05)); 
epsilon = 0.1;

N = 100000;
x = zeros(1,N);
J = zeros(1,N); 

for i = 1:N
    
 J_atual = x_atual ^ 2;
 x(i) = x_atual;
 J(i) = J_atual;
 
 r = random('unif', -1, 1);
 x_futuro = x_atual + epsilon * r;
 J_futuro = x_futuro ^ 2;
 J_dif = J_futuro - J_atual;

 if J_dif < 0
     x_atual = x_futuro;
 else
     beta = rand(1,1);

     if beta < exp(-(J_dif)/0.1)
         x_atual = x_futuro;
     end
 end

end

% Plot dos gráficos

hist(x((N-10000):end), 100)
set(gca, 'FontSize', 30)
title('Histograma dos últimos 10000 estados x(n)', 'FontSize', 30)
disp(['Média dos últimos 10000 estados x(n): ' num2str(mean(x((N-10000):end)))])
disp(['Variância dos últimos 10000 estados x(n): ' num2str(var(x((N-10000):end)))])

% Valores encontrados para média e variância
% 
% Média dos últimos 10000 estados x(n): -0.014273
% Variância dos últimos 10000 estados x(n): 0.050015

figure

hist(J((N-10000):end), 100)
set(gca, 'FontSize', 30)
title('Histograma dos últimos 10000 valores J(x)', 'FontSize', 30)

figure

plot(x, 'o ')
set(gca, 'FontSize', 30)
xlabel('n','FontSize', 30) 
ylabel('x(n)', 'FontSize', 30)
title('Valores dos estados x(n)', 'FontSize', 30) 

figure

plot(1:10000, x(1:10000), 'o ')
set(gca, 'FontSize', 30)
xlabel('n', 'FontSize', 30) 
ylabel('x(n)', 'FontSize', 30) 
title('Valores dos primeiros 10000 estados x(n)', 'FontSize', 30) 
axis([1 10000 -0.8 0.8])

figure 

plot(90000:100000, x((N-10000):end), 'o ') 
xlabel('n', 'FontSize', 30) 
ylabel('x(n)', 'FontSize', 30) 
title('Valores dos últimos 10000 estados x(n)', 'FontSize', 30) 
axis([90000 100000 -0.8 0.8]) 