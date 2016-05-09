%% Montagem da matriz de transição M

C1 = [0 1/4 1/4 1/4 1/4]';
C2 = [exp(-3)/4 (3-exp(-3)-exp(-2)-exp(-1))/4 exp(-1)/4 1/4 exp(-2)/4]';
C3 = [exp(-2)/4 1/4 (2-exp(-2)-exp(-1))/4 1/4 exp(-1)/4]';
C4 = [exp(-4)/4 exp(-1)/4 exp(-2)/4 (4-exp(-4)-exp(-3)-exp(-2)-exp(-1))/4 exp(-3)/4]';
C5 = [exp(-1)/4 1/4 1/4 1/4 (1-exp(-1))/4]';

M = [C1 C2 C3 C4 C5];

%% Inicialização p0

p0 = zeros(1,5)';

% Geração aleatória de p0

r = rand();
b = 1;

for i = 1:(length(p0)-1)
    p0(i) = r;
    c = b - r;
    p0(i+1) = c;
    r = random('unif', 0, c);
    b = c;
end

%% Cálculo do vetor invariante

p_atual = p0;
contador = 1;

while(true)
    p_seguinte = M * p_atual;
    if norm((p_seguinte - p_atual), 2) == 0
        break
    end
    p_atual = p_seguinte;
    contador = contador + 1;
end

p_inv = p_atual
contador

% Resultados
% 
% p_inv =
% 
%     0.0117
%     0.2341
%     0.0861
%     0.6364
%     0.0317
% 
% 
% contador =
% 
%     84