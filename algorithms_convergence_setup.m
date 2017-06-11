% Runs algorithms convergence. No loops to speed up parfor computattions

%% Constant parameters
eps = 0.01;
num_trials = 50;
num_iter = 100;
precision = 1e-6;
write_precision = 6;
graph_size = 50;
prob = 0.5;

folder = 'convergence_results/';
mkdir(folder);

%% Choose methods
do_irls = true;
do_grad = true;
do_singval = true;
do_logdet = true;

%% Run scripts

disp('Running p=0.5 and q=0.8');
p = 0.5;
q = 0.8;

algorithms_convergence

do_logdet = false;

disp('Running p=0.1 and q=0.5 (no logdet');
p = 0.1;
q = 0.5;

algorithms_convergence
