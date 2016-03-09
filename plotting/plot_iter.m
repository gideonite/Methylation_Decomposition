filename = 'hmm1.0_jerry1';
load(['results/' filename '.mat']);
h = figure;
x = 1:length(curve_errors);
plot(x, curve_errors, 'b')
xlabel('iterations')
ylabel('MSE')
% legend('best_ce','decomp_errors', 'Location','northeast')
