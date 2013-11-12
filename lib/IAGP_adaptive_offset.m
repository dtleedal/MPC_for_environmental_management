function [g_k, p_k] = IAGP_adaptive_offset( g_k_1, p_k_1, y_hat, y, nvr)

% p_k_k_1 = p_k_1 + nvr;
% p_k = p_k_k_1 - (p_k_k_1.^2 * y_hat.^2)./(1 + p_k_k_1 * y_hat.^2);
% %g_k = g_k_1 + p_k * y_hat * (y - g_k_1 * y_hat); % the adaptive gain
% g_k = g_k_1 + p_k * y_hat * (y - (g_k_1 * y_hat));
% offset = (g_k*y_hat) - y_hat;


p_k_k_1 = p_k_1 + nvr;
p_k = p_k_k_1 - (p_k_k_1.^2 * y_hat.^2)./(1 + p_k_k_1 * y_hat.^2);
g_k = g_k_1 + p_k * y_hat * (y - (g_k_1 + y_hat));
if g_k < -2
    warning('Adaptive gain less than 0.5. Setting to 0.5')
    g_k = -2;
end
if g_k > 2
    warning('Adaptive gain greater than 1.5. Setting to 1.5')
    g_k = 2;
end
