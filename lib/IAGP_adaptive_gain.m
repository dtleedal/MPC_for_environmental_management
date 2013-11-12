function [g_k, p_k] = IAGP_adaptive_gain( g_k_1, p_k_1, y_hat, y, nvr)

p_k_k_1 = p_k_1 + nvr;
p_k = p_k_k_1 - (p_k_k_1.^2 * y_hat.^2)./(1 + p_k_k_1 * y_hat.^2);
g_k = g_k_1 + p_k * y_hat * (y - g_k_1 * y_hat);
if g_k < 0.1
    warning('Adaptive gain less than 0.5. Setting to 0.5')
    g_k = 0.1;
end
if g_k > 1.9
    warning('Adaptive gain greater than 1.5. Setting to 1.5')
    g_k = 1.9;
end

end

