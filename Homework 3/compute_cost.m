%% Che Jin Goh | UID: 905724540

function cost = compute_cost(AL, Y)
    % Computes the cross-entropy loss
    % AL: Predicted outputs (probabilities), shape (number of classes, number of examples)
    % Y: True labels, shape (number of classes, number of examples)
    % Returns: Scalar cost (average cross-entropy loss)

    m = size(Y, 2); % Number of examples
    cost = -sum(sum(Y .* log(AL + 1e-8))) / m; % Add epsilon to prevent log(0)
end
