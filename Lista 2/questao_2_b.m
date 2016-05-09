%% Montagem da matriz de transição M

C1 = [0 1/4 1/4 1/4 1/4]';
C2 = [exp(-3)/4 (3-exp(-3)-exp(-2)-exp(-1))/4 exp(-1)/4 1/4 exp(-2)/4]';
C3 = [exp(-2)/4 1/4 (2-exp(-2)-exp(-1))/4 1/4 exp(-1)/4]';
C4 = [exp(-4)/4 exp(-1)/4 exp(-2)/4 (4-exp(-4)-exp(-3)-exp(-2)-exp(-1))/4 exp(-3)/4]';
C5 = [exp(-1)/4 1/4 1/4 1/4 (1-exp(-1))/4]';

M = [C1 C2 C3 C4 C5];

%% Cálculo dos estados

X = zeros(1,4);
T = 4;

p0 = [1 0 0 0 0]';

p = p0;

for i = 0:(T-1)
    if i==0 % Estado inicial
        X(i+1) = 1;
    else
        p = M * p;
        r = rand();
        p_cdf = cumsum(p);
        
        for j = 1:length(p_cdf)
            if r<=p_cdf(j)
                switch j
                    case 1
                        x = 1;
                    case 2
                        x = 2;
                    case 3
                        x = 3;
                    case 4
                        x = 4;
                    case 5
                        x = 5;
                end
                break
            end
        end
        
        X(i+1) = x;
    end
end

X