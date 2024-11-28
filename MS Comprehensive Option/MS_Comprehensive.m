clear; close all; clc;

% Generate Data
rng(42); % Seed for reproducibility
x_generated = linspace(0, 5, 10);
n_true = 0.06;
a_true = 0.25;
m_true = 0.57;
b_true = 0.11;
noise = 0.001 * randn(size(x_generated));
y_generated = n_true * exp(-a_true * (m_true * x_generated + b_true).^2) + noise;

% Define Learning Rates and Epochs for Experiments
learning_rates = [0.0001, 0.001, 0.01]; % Learning rates to test
N = length(x_generated); % Number of data points

% Initialize storage for loss histories and final parameters
loss_histories_lr = cell(length(learning_rates), 1);
final_parameters_lr = zeros(length(learning_rates), 2); % [m, b]

%% Vary Learning Rates
for lr_idx = 1:length(learning_rates)
    learning_rate = learning_rates(lr_idx);

    % Initialize Parameters
    m = 0; 
    b = 0; 
    epochs = 10000; 
    loss_history = zeros(epochs, 1);

    % Gradient Descent
    for epoch = 1:epochs
        % Compute predictions
        y_pred = m * x_generated + b;

        % Compute gradients
        dL_dm = (-2 / N) * sum(x_generated .* (y_generated - y_pred));
        dL_db = (-2 / N) * sum(y_generated - y_pred);

        % Update parameters
        m = m - learning_rate * dL_dm;
        b = b - learning_rate * dL_db;

        % Compute loss (MSE)
        loss = mean((y_generated - y_pred).^2);
        loss_history(epoch) = loss;
    end

    % Store results
    loss_histories_lr{lr_idx} = loss_history;
    final_parameters_lr(lr_idx, :) = [m, b];
end

%% Vary Epochs 
learning_rate = 0.001;

% Initialize Parameters
m = 0; % Initial guess for slope
b = 0; % Initial guess for intercept
loss_history = zeros(epochs, 1);

% Gradient Descent
for epoch = 1:epochs
    % Compute predictions
    y_pred = m * x_generated + b;

    % Compute gradients
    dL_dm = (-2 / N) * sum(x_generated .* (y_generated - y_pred));
    dL_db = (-2 / N) * sum(y_generated - y_pred);

    % Update parameters
    m = m - learning_rate * dL_dm;
    b = b - learning_rate * dL_db;

    % Compute loss (MSE)
    loss = mean((y_generated - y_pred).^2);
    loss_history(epoch) = loss;
end

% Store results
loss_history_epochs = loss_history;
final_parameters_epochs = [m, b];

%% Plot Results

% Plot 1: Fitted Line vs. Actual Data
figure()
scatter(x_generated, y_generated, 'b', 'filled', 'DisplayName', 'Actual Data');
hold on
m = final_parameters_lr(2, 1);
b = final_parameters_lr(2, 2);
plot(x_generated, m * x_generated + b, 'r', 'LineWidth', 1.5, 'DisplayName', 'Fitted Line');
title('Fitted Line vs. Actual Data')
xlabel('x')
ylabel('y')
legend('show')
grid on

% Plot 2: Loss vs. Epochs for Different Learning Rates
figure()
hold on
for lr_idx = 1:length(learning_rates)
    plot(loss_histories_lr{lr_idx}, 'LineWidth', 1.5,'DisplayName', sprintf('LR = %.4f', learning_rates(lr_idx)));
end
title('Loss vs. Epochs for Different Learning Rates')
xlabel('Epochs')
ylabel('MSE Loss')
legend('show')
grid on

% Plot 3: Loss vs. Epochs for the Largest Epoch Configuration
figure()
plot(1:epochs, loss_history_epochs, 'r', 'LineWidth', 1.5, 'DisplayName', '10,000 Epochs');
title('Loss vs. Epochs for 10,000 Epochs');
xlabel('Epochs')
ylabel('MSE Loss')
legend('show')
grid on

% Display Final Parameters
fprintf('Final Parameters for Learning Rate Experiment:\n');
for lr_idx = 1:length(learning_rates)
    fprintf('LR = %.4f: m = %.6f, b = %.6f\n', learning_rates(lr_idx), ...
        final_parameters_lr(lr_idx, 1), final_parameters_lr(lr_idx, 2));
end
fprintf('\nFinal Parameters for Epoch Experiment:\n');
fprintf('m = %.6f, b = %.6f\n', final_parameters_epochs(1), final_parameters_epochs(2));
