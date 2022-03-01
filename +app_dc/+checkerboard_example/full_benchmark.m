% This script runs full DC checkerboards benchmark
% to stress the Woodbury-MINRES-AMG solver. This takes
% too much time and memory and such is not tested on CI.
%
% Postprocess by other script in this directory as
% appropriate.

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2ddir-n0004', 'direct');  % ensure warm start
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2ddir-n0004', 'direct');
timings(2) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2ddir-n0008', 'direct');
timings(3) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2ddir-n0016', 'direct');
timings(4) = app_dc.checkerboard_example.solve(2,   32, 1, 2, '2ddir-n0032', 'direct');
%timings(5) = app_dc.checkerboard_example.solve(2,   64, 1, 2, '2ddir-n0064', 'direct');
%timings(6) = app_dc.checkerboard_example.solve(2,  128, 1, 2, '2ddir-n0128', 'direct');
%timings(7) = app_dc.checkerboard_example.solve(2,  256, 1, 2, '2ddir-n0256', 'direct');
save('checkerboard-timings-2ddir.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dfw-n0004', 'krylov-full-woodbury');  % ensure warm start
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dfw-n0004', 'krylov-full-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dfw-n0008', 'krylov-full-woodbury');
timings(3) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dfw-n0016', 'krylov-full-woodbury');
timings(4) = app_dc.checkerboard_example.solve(2,   32, 1, 2, '2dfw-n0032', 'krylov-full-woodbury');
timings(5) = app_dc.checkerboard_example.solve(2,   64, 1, 2, '2dfw-n0064', 'krylov-full-woodbury');
%timings(6) = app_dc.checkerboard_example.solve(2,  128, 1, 2, '2dfw-n0128', 'krylov-full-woodbury');
%timings(7) = app_dc.checkerboard_example.solve(2,  256, 1, 2, '2dfw-n0256', 'krylov-full-woodbury');
%timings(8) = app_dc.checkerboard_example.solve(2,  512, 1, 2, '2dfw-n0512', 'krylov-full-woodbury');
%timings(9) = app_dc.checkerboard_example.solve(2, 1024, 1, 2, '2dfw-n1024', 'krylov-full-woodbury');
save('checkerboard-timings-2dfw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dnw-n0004', 'krylov-no-woodbury');  % ensure warm start
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dnw-n0004', 'krylov-no-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dnw-n0008', 'krylov-no-woodbury');
%timings(3) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dnw-n0016', 'krylov-no-woodbury');
%timings(4) = app_dc.checkerboard_example.solve(2,   32, 1, 2, '2dnw-n0032', 'krylov-no-woodbury');
%timings(5) = app_dc.checkerboard_example.solve(2,   64, 1, 2, '2dnw-n0064', 'krylov-no-woodbury');
save('checkerboard-timings-2dnw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(3, [2, 2], 2, 2, '3ddir-2x2', 'direct');
%timings(2) = app_dc.checkerboard_example.solve(3, [2, 3], 2, 2, '3ddir-2x3', 'direct');
%timings(3) = app_dc.checkerboard_example.solve(3, [2, 4], 2, 2, '3ddir-2x4', 'direct');
save('checkerboard-timings-3ddir.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(3, [2, 2], 2, 2, '3dfw-2x2', 'krylov-full-woodbury');
timings(2) = app_dc.checkerboard_example.solve(3, [2, 3], 2, 2, '3dfw-2x3', 'krylov-full-woodbury');
%timings(3) = app_dc.checkerboard_example.solve(3, [2, 4], 2, 2, '3dfw-2x4', 'krylov-full-woodbury');
%timings(4) = app_dc.checkerboard_example.solve(3, [2, 5], 2, 2, '3dfw-2x5', 'krylov-full-woodbury');
save('checkerboard-timings-3dfw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(3, [2, 2], 2, 2, '3dnw-2x2', 'krylov-no-woodbury');
%timings(2) = app_dc.checkerboard_example.solve(3, [2, 3], 2, 2, '3dnw-2x3', 'krylov-no-woodbury');
%timings(3) = app_dc.checkerboard_example.solve(3, [2, 4], 2, 2, '3dnw-2x4', 'krylov-no-woodbury');
save('checkerboard-timings-3dnw.mat', 'timings');
