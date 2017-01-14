sample_minmax = [0  0   5   3.4   5.6
        4  4  34   4.6   52  ];
s = SampleNdRange(sample_minmax)
s.sample(10)

sample_step = [1 0 3 1.2 0];   %Zero implies continous sample
s2 = SampleNdRange(sample_minmax, sample_step)
s2.sample(10)

s3 = SampleNdRange(sample_minmax, sample_step, 'sobol')  %Default is psuedo-random
s3.sample(10) %10 randomly scrambled 5-D sobol samples over range
