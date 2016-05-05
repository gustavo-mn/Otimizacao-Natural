X = (1:5)';
J = [0.5 0.2 0.3 0.1 0.4]';

T = 0.1;

for i = 1:length(J)
    B(i) = exp(-(J(i)/T));
end

B = B';

PB = B/sum(B);

tabela = [X J B PB]