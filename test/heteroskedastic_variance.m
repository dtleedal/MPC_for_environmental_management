function het_var = heteroskedastic_variance( y,n )
% n = 9;
% y = randn(40,1)
% y is the input series
% n is the window width
het_var = nan(size(y));
for i = n+1:length(y)-n
    het_var(i) = sum((y(i-n:i+n).^2));%std(y(i-n:i+n));
end


end

