function plot_dados_questao_2(N_ger, abordagem, memetico)

executa = 1;

if memetico

    if exist(['dados_questao_2_' num2str(N_ger) '_geracoes_' abordagem '.mat'], 'file')
        load(['dados_questao_2_' num2str(N_ger) '_geracoes_' abordagem '.mat'])
    else
        executa = 0;
    end

else
    
    if exist('dados_questao_2_nao_memetico.mat', 'file')
        load('dados_questao_2_nao_memetico.mat')
    else
        executa = 0;
    end
    
end

N = zeros(1,length(geracoes_otimas));

if executa
    
    for t = 1:length(geracoes_otimas)       
        N(t) = geracoes_otimas(t).n;
        
        clc
        
        disp(['Execu��o: t = ' num2str(t)])
        
        figure(1)
        
        stem(geracoes_otimas(t).maximos);
        axis([0 length(geracoes_otimas(t).maximos) 0 25])
        xlabel('Gera��o', 'FontSize', 30);
        ylabel('Fitness M�xima', 'FontSize', 30);
        title('Melhores fitness por gera��o', 'FontSize', 30);
        set(gca, 'FontSize', 30);
        
        figure(2)

        stem(geracoes_otimas(t).minimos);
        axis([0 length(geracoes_otimas(t).minimos) 0 25])
        xlabel('Gera��o', 'FontSize', 30);
        ylabel('Fitness M�nima', 'FontSize', 30);
        title('Piores fitness por gera��o', 'FontSize', 30);
        set(gca, 'FontSize', 30);
        
        figure(3)

        stem(geracoes_otimas(t).medias);
        axis([0 length(geracoes_otimas(t).medias) 0 25])
        xlabel('Gera��o', 'FontSize', 30);
        ylabel('M�dia Fitness', 'FontSize', 30);
        title('M�dia dos fitness por gera��o', 'FontSize', 30);
        set(gca, 'FontSize', 30);
        
        pause();
    end
    
    disp(['M�dia de gera��es para alcan�ar o m�ximo: ' num2str(mean(N))])
    disp(['Desvio-padr�o de gera��es para alcan�ar o m�ximo: ' num2str(std(N, 1))])
else
    disp('N�o h� dados dispon�veis')
end

end