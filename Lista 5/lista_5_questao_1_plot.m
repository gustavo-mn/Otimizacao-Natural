function lista_5_questao_1_plot()

clear all; close all; clc

load('dados_lista_5_questao_1.mat')

stem(melhor_fitness_execucao)
xlabel('Execução', 'FontSize', 30)
ylabel('Melhor fitness', 'FontSize', 30)
title('Melhor fitness por execução', 'FontSize', 30)
set(gca, 'FontSize', 30)

figure;
 
% edges = [1e-5 1e-3 0.1 1 1.5 2 2.5 3];
% N = histc(melhor_fitness_execucao, edges);
% bar(edges, N, 'histc')
hist(melhor_fitness_execucao)
xlabel('Melhor fitness', 'FontSize', 30)
title('Histograma das melhores fitness', 'FontSize', 30)
set(gca, 'FontSize', 30)

figure;

stem(geracao_otima_execucao)
xlabel('Execução', 'FontSize', 30)
ylabel('Geração', 'FontSize', 30)
title('Gerações ótimas por execução', 'FontSize', 30)
set(gca, 'FontSize', 30)

figure;

hist(geracao_otima_execucao)
xlabel('Geração', 'FontSize', 30)
title('Histograma das melhores gerações', 'FontSize', 30)
set(gca, 'FontSize', 30)
