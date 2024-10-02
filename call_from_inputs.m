inputs=load("C:\Users\gwosiak\Documents\MATBOX tests\Size distribution\Inputs.mat");


%inputs.phase.volumefraction.along_3rd_axis;

%vartest=inputs.phase.size_histogram.along_3rd_axis; %% not working 

inputs.phase.size_histogram.along_3rd_axis=vartest;

run_number=1;

savefolder='C:\Users\gwosiak\Documents\MATBOX tests\Size distribution\';

  save_progression = false; 
  save_verification = true;

Microstructure_generation_stochastic_frominputs(inputs, 1, run_number,save_progression, save_verification);

