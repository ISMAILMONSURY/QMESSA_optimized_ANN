net = newff(i,t,[10],{'tansig','purelin'},'trainlm','learngdm','mse');
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
net_lm = train(net,i,t);