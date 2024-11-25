%% Che Jin Goh | UID: 905724540

function Z = linear_forward(A_prev, W, b)
% Performs the linear part of the forward propagation for a single layer
% A_prev: Activations from the previous layer, shape (size of previous layer, number of examples)
% W: Weight matrix of the current layer, shape (size of current layer, size of previous layer)
% b: Bias vector of the current layer, shape (size of current layer, 1)
% The method computes Z = W * A_prev + b, where Z is the input to the activation function.
% Returns: Z, the pre-activation values, shape (size of current layer, number of examples)

    Z = W * A_prev + b;
end