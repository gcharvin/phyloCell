%%

outputDirectory = 'png20120823';

%%

convertImagesToPng8('L:\common\movies\steffen\2011\201101\20110121', 'bud9_delta_mut-1', 1, 1, 6, outputDirectory);

%%

base = 'L:\common\movies\steffen\2012\201207\120726_CP03_deflector_regular_setup';
project = '120726_CP03_deflector';

for i=5:29
    fprintf('%d/%d\n', i, 29);
    tic
    convertImagesToPng8(base, project, i, 1, 6, outputDirectory);
    convertImagesToPng8(base, project, i, 3, 0, outputDirectory);
    toc
end

%

base = 'L:\common\movies\steffen\2012\201204\120405_CP03_full_spectum';
project = '120405_CP03_full_spectrum_hi-pH';
n = 28;

for i=3:n
    fprintf('%d/%d\n', i, n);
    tic
    convertImagesToPng8(base, project, i, 1, 6, outputDirectory);
    convertImagesToPng8(base, project, i, 3, 0, outputDirectory);
    toc
end

%

base = 'L:\common\movies\steffen\2012\201207\120720_deflector_3min_intervals';
project = '120720_deflector_3min_intervals';
n = 10;

for i=3:n
    fprintf('%d/%d\n', i, n);
    tic
    convertImagesToPng8(base, project, i, 1, 6, outputDirectory);
    convertImagesToPng8(base, project, i, 2, 0, outputDirectory);
    toc
end

%%

base = 'L:\common\movies\youlian\2011\nov\25112011';
project = 'SRX1-GFP_pedigree_2';
n = 0;%4;

for i=1:n
    fprintf('%d/%d\n', i, n);
    tic
    convertImagesToPng8(base, project, i, 1, 6, outputDirectory);
    toc
end

%%

base = 'L:\common\movies\steffen\2011\201109\110930_CP03-1';
project = 'Stress_and_cell_cycle';

for i=[2 5 31]
    fprintf('%d/%d\n', i, 31);
    tic
    convertImagesToPng8(base, project, i, 1, 6, outputDirectory);
    convertImagesToPng8(base, project, i, 3, 0, outputDirectory);
    toc
end
