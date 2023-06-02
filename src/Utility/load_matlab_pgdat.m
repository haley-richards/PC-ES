function [pgdat]=load_matlab_pgdat(valpth)
basepth=[valpth];
fname='coordinates.mat';
pgdat=load([loadpth, fname]);
end