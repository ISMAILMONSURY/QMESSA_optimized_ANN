% main_QMESSABP.m
% Code generated by MD. Ismail Monsury
% How to implement optimization algorithm(QMESSA) in neural network
% Main purpose is to optimize weight & bias values to finetune the model
% and avoid local minima problem.
% INITIALIZE THE NEURAL NETWORK  %

% inputs for the neural net
inputs = i;

% targets for the neural net
targets = t;

% number of neurons
n = 10;

% create a neural network
net = feedforwardnet(n);

% configure the neural network for this dataset
net = configure(net, inputs, targets);

% get the normal NN weights and bias
getwb(net);

% error MSE normal NN
error = targets - net(inputs);
calc = mean(error.^2)/mean(var(targets',1));

% create handle to the MSE_TEST function, that
% calculates MSE

h = @(bestX) NMSE(bestX, net, inputs, targets);



% running the improved sparrow search  optimization algorithm with desired options
[bestX, err_SSA] = QMESSA(14*n+n+1,h )
net = setwb(net, bestX');

% get the SSA optimized NN weights and bias(mahi)
getwb(net);

% error MSE SSA optimized NN
error = targets - net(inputs);
calc = mean(error.^2)/mean(var(targets',1));