%% Che Jin Goh | UID: 905724540

%% Main File
clear all; close all; clc;

%% Load Training and Testing Data
[X_train, Y_train, X_test, Y_test] = load_train_and_test_data();

%% Define fixed hyperparameters
input_size = size(X_train, 1); % Input size of the FNN
output_size = size(Y_train, 1); % Output size (number of classes)
neurons = 64; % Neurons in each layer
epochs = 150; % Fixed number of epochs
lr = 0.01; % Fixed learning rate
numLayerConfigs = [2, 3, 5]; % Number of layers to test

%% Initialize storage for combined plots
allTrainLoss = cell(length(numLayerConfigs), 1);
allTestAccuracy = cell(length(numLayerConfigs), 1);

%% Loop through each number of layer configuration
for layerIndex = 1:length(numLayerConfigs)
    numLayer = numLayerConfigs(layerIndex);
    fprintf('Running model with %d layers...\n', numLayer);
    
    % Define layer dimensions
    layer_dims = zeros(1, numLayer + 2);
    layer_dims(1) = input_size;
    layer_dims(end) = output_size;
    for i = 1:numLayer
        layer_dims(i+1) = neurons;
    end
    
    % Initialize parameters
    parameters = initialize_parameters(layer_dims);
    
    % Initialize progress trackers
    m = size(X_train, 2);
    batch_size = 64;
    num_batches = floor(m / batch_size);
    trainLoss = zeros(epochs, 1);
    testAccuracy = zeros(epochs, 1);

    % Training loop
    for epoch = 1:epochs
        % Shuffle the training data
        indices = randperm(m);
        X_train = X_train(:, indices);
        Y_train = Y_train(:, indices);
        
        % Train in mini-batches
        cost = zeros(num_batches, 1);
        for j = 1:num_batches
            X_batch = X_train(:, (j-1)*batch_size+1:j*batch_size);
            Y_batch = Y_train(:, (j-1)*batch_size+1:j*batch_size);
            forward_pass = forward_propagation(X_batch, parameters);
            cost(j) = compute_cost(forward_pass{end}, Y_batch); % Compute scalar cost
            gradients = backward_propagation(X_batch, Y_batch, parameters, forward_pass);
            parameters = update_parameters(parameters, gradients, lr);
        end
        
        % Evaluate on the test set
        Y_pred = predict(X_test, parameters);
        acc = accuracy(Y_pred, Y_test);
        
        % Log training loss and testing accuracy
        trainLoss(epoch) = mean(cost);
        testAccuracy(epoch) = acc;
        
        % Print progress
        fprintf('Epoch %d/%d: Loss = %f, Accuracy = %f\n', epoch, epochs, trainLoss(epoch), acc);
    end
    
    %% Store results for combined plot
    allTrainLoss{layerIndex} = trainLoss;
    allTestAccuracy{layerIndex} = testAccuracy;
end

%% Combined Plot for Training Loss and Testing Accuracy
figure;

% Plot training loss
subplot(2, 1, 1);
hold on;
for layerIndex = 1:length(numLayerConfigs)
    numLayer = numLayerConfigs(layerIndex);
    plot(1:epochs, allTrainLoss{layerIndex}, 'DisplayName', sprintf('%d Layers', numLayer));
end
title('Training Loss vs. Epochs');
xlabel('Epochs');
ylabel('Training Loss');
legend('show');
grid on;

% Plot testing accuracy
subplot(2, 1, 2);
hold on;
for layerIndex = 1:length(numLayerConfigs)
    numLayer = numLayerConfigs(layerIndex);
    plot(1:epochs, allTestAccuracy{layerIndex}, 'DisplayName', sprintf('%d Layers', numLayer));
end
title('Testing Accuracy vs. Epochs');
xlabel('Epochs');
ylabel('Testing Accuracy');
legend('show');
grid on;

% Automatically save the combined figure
saveas(gcf, 'task3_combined_num_layers.png');
fprintf('Figure saved as task3_combined_num_layers.png\n');
