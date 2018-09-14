clear; clc;

load mydata_180908;
x = double(data.x1'); y = double(data.y1');
data_size = 6000;
train_ratio = 0.2;
% train_x = double(data.train_x'); train_y = double(data.train_y');
% test_x = double(data.test_x'); test_y = double(data.test_y');
index = 5000;
train_x = x(:,1:index); train_y = y(:,1:index);
test_x = x(:,index+1:6000); test_y = y(:,index+1:6000);

% rand_x_lost = randn(60,2000)-0.5;
% train_x = train_x - 0.1*rand_x_lost;

hyperParas.debug = 0;

hyperParas.arch = [32, 256, 128, 64, 32, 24];
hyperParas.numLayer = numel(hyperParas.arch);
hyperParas.outDim = hyperParas.arch(end);
hyperParas.actFunc = 'relu'; % sigm, tanh, relu, 1
hyperParas.loss = 'square'; % crossEnt, square

hyperParas.learnRate = 0.075;
hyperParas.batchSize = 200; % trainSize 必须是 batchSize 的整数倍，否则最后一个batch会出错
hyperParas.numEpochs = 300;

modelParas = nninit(hyperParas);
[modelParas, losses] = nntrain(hyperParas, modelParas, train_x, train_y);
err = nntest(hyperParas, modelParas, test_x, test_y);
figure = figure('color',[1,1,1]); plot(losses(1:1:end));