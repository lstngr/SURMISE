function out=read_output( filename )
% READ_OUPUT reads a SURMISE output file and stores its contents in an
% appropriate data structure.
%
% Arguments:
%   - filename: String with a path to the SURMISE output
%
% Returns:
%   - out: Structure storing the simulation's output

% === NOTE ===
% Will be adapted to CSV formatted files in a future version. Will expect
% following structure (header part, data part)
% " property1,value1
% " property2,value2
% " ...
% " propertyN,valueN
% " 0,x,y,vx,vy,fx,fy
% " 1,x,y,vx,vy,fx,fy
% " ...
% " M,x,y,vx,vy,fx,fy
% Where M,N known integers (particle count and parameter number)
% === ==== ===

allfiles = dir([filename,'.leafs.*']);

data = csvread([allfiles(1).folder,'/',allfiles(1).name]);
pnum = max(data(:,1)) + 1;
disp([num2str(pnum),' particles found'])
out = zeros(length(allfiles),pnum,7);
% Indexes of out refer to
% - time step
% - particle index
% - variable of interest (mass, pos_x/y, vel_x/y, frc_x/y)

for ifile=1:length(allfiles)
    data = csvread([allfiles(ifile).folder,'/',allfiles(ifile).name]);
    for il=1:size(data,1)
        out(ifile,data(il,1)+1,1:7) = data(il,2:8);
%         out(ifile,data(il,1)+1,2) = data(il,3);
%         out(ifile,data(il,1)+1,3) = data(il,4);
%         out(ifile,data(il,1)+1,4) = data(il,5);
%         out(ifile,data(il,1)+1,5) = data(il,6);
%         out(ifile,data(il,1)+1,6) = data(il,7);
%         out(ifile,data(il,1)+1,7) = data(il,8);
    end
end
end