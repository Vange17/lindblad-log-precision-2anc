function verify_variance(output_dir)
% Small regression checks; output_dir must not contain previous test outputs.
variance_dir=fileparts(mfilename('fullpath'));
assert(isfolder(output_dir),'Create the test output directory first.');

% The first draw with this seed enters kA>0. The upstream initialization
% fails with undefined xlab_prod; the corrected sampler must complete.
load(fullfile(fileparts(variance_dir),...
 'miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat'),'mundt1');
rng(1989,'twister');
assert(rand>mundt1(1)/sum(mundt1));
run_variance(2,1989,fullfile(output_dir,'regression_fixed'));

% The fixture was evaluated using expm and the original ancilla dilation
% at every step, with the same approved initialization correction.
stem=fullfile(output_dir,'equivalence_optimized');
run_variance(20,20260907,stem);
actual=load([stem '.mat']);
reference=load(fullfile(variance_dir,'tests','reference_20.mat'));
assert(actual.num_NMt==reference.num_NMt && actual.seed==reference.seed);
assert(isequal(actual.t_vec,reference.t_vec));
trajectory_max_difference=max(abs(actual.values-reference.values),[],'all');
assert(trajectory_max_difference<1e-11);
assert(actual.baseline_max_difference<1e-11);
assert(all(actual.values(:,1)==1));
assert(all(isfinite(actual.values),'all'));
assert(all(actual.trajectory_variance>=0));
assert(max(abs(actual.trajectory_std.^2-actual.trajectory_variance))<1e-12);
assert(max(abs(actual.standard_error_N5million.^2*5e6-actual.trajectory_variance))<1e-12);

% Frozen files must not be overwritten by a repeated invocation.
refused_overwrite=false;
try
 run_variance(20,20260907,stem);
catch exception
 assert(contains(exception.message,'Output already exists'),exception.message);
 refused_overwrite=true;
end
assert(refused_overwrite);
fprintf('PASS: initialization, reference trajectories, variance/SEM, overwrite protection.\n');
fprintf('Maximum trajectory difference: %.4g\n',trajectory_max_difference);
end
