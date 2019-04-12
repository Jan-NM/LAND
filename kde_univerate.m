x = [-3.9 -2.9 -1.8 0.9 1.5 2];
histogram(x,8)
hold on
plot(x,0,'o')
hold off
% define kernel bandwith
h = 0.75;
% sampling rate
samplingRate = 100;
% generate line vector
y = linspace(-5, 3, samplingRate);
% kernel function
K = 0;
for i=1:numel(x)
    K = 1/h * 1/sqrt(2*pi) .* exp(-0.5.*((y-x(i))./h).^2) + K; 
end

K = 1/numel(x) .* K;

plot(y, K);
% test single values
hold on
K = 0;
for i=1:numel(x(1))
    K = 1/h * 1/sqrt(2*pi) .* exp(-0.5.*((y-x(i))./h).^2) + K;
end
% K = 1/numel(x(1)) .* K;
K = 1/numel(x) .* K;
plot(y, K);