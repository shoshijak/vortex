clear
close all
clc

filename = 'diegoBinaryN12000';
fid = fopen(filename, 'r');

nSources = fread(fid,1,'integer*4')
xSources = fread(fid, nSources, 'real*8');
ySources = fread(fid, nSources, 'real*8');
sources = fread(fid, nSources, 'real*8');

nTargets = fread(fid,1,'integer*4')
xTargets = fread(fid, nTargets, 'real*8');
yTargets = fread(fid, nTargets, 'real*8');
pressure = fread(fid, nTargets, 'real*8');

if(~isempty(fread(fid)))
    error('Complete file not read in')
end
