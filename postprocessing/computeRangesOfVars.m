infile = '../model/output/test_model2_long.mat';
load(infile)

transient = 10000;

Y = Y(:,transient+1:end);

num_eqns_per_vertex = 5;
num_eqns_per_edge = 1;
neuron_classes = 4;

neuron_mean = zeros(num_eqns_per_vertex, neuron_classes);
neuron_std = zeros(num_eqns_per_vertex, neuron_classes);
neuron_min = zeros(num_eqns_per_vertex, neuron_classes);
neuron_max = zeros(num_eqns_per_vertex, neuron_classes);
edge_mean = zeros(1, neuron_classes);

for i = 1:neuron_classes
    %% vertex variables
    firstvar=(i-1)*num_eqns_per_vertex+1;
    idx = firstvar:firstvar+num_eqns_per_vertex-1;
    disp(idx')
    neuron_mean(:,i) = mean(Y(idx, :), 2);
    neuron_std(:,i)  = std(Y(idx, :), 0, 2);
    neuron_min(:,i) = min(Y(idx,:), [], 2);
    neuron_max(:,i) = max(Y(idx,:), [], 2);
end

disp( min( neuron_min, [], 2 ) );
disp( max( neuron_max, [], 2 ) );

    