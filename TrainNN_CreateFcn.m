%Nathan Parrish
%CCAR 2017
%
%This function is designed to be called from Julia. We pass in the training
%data and get a trained network as a result, with the linear scale/shift
%terms output from the function.
%

function scaleShift = TrainNN_CreateFcn(NN_input, NN_output, NN_size, fcn_name)

%fix the format of NN_size (the training function is very picky about this)
NN_size = NN_size(:)'; %ensure that NN_size is a row vector
NN_size = double(NN_size); %ensure that NN_size is a double

%% Neural Network training
M_train = length(NN_output);

%scale on inputs:
scale_in = std(NN_input,0,2);
shift_in = mean(NN_input, 2);

%scale on outputs:
scale_out = std(NN_output,0,2);
shift_out = mean(NN_output,2);

x = (NN_input - repmat(shift_in, 1, M_train)) ./ repmat(scale_in, 1, M_train);
t = (NN_output - repmat(shift_out, 1, M_train)) ./ repmat(scale_out, 1, M_train);

net = fitnet(NN_size);

net.trainParam.max_fail = 50;
net.trainParam.epochs = 2000;
% net.divideFcn = 'divideblock';
% net.trainFcn = 'trainscg';

[net,tr] = train(net,x,t,'useParallel','yes');
% [net,tr] = train(net,x,t,'useGPU','only');

%Generate function that we can call from Julia:
genFunction(net, fcn_name);

if length(scale_in) > length(scale_out) %bigger inputs than outputs
    %pad with NaN's
    scale_out = [scale_out; NaN * ones(length(scale_in) - length(scale_out), 1)];
    shift_out = [shift_out; NaN * ones(length(shift_in) - length(shift_out), 1)];

elseif length(scale_in) < length(scale_out) %bigger outputs than inputs
    %pad with NaN's
    scale_in = [scale_in; NaN * ones(length(scale_out) - length(scale_in), 1)];
    shift_in = [shift_in; NaN * ones(length(shift_out) - length(shift_in), 1)];

end

scaleShift = [scale_in, shift_in, scale_out, shift_out];

end
