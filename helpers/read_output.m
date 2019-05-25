function out=read_output( filename )
% READ_OUPUT reads a SURMISE output file and stores its contents in an
% appropriate data structure.
%
% Arguments:
%   - filename: String with a path to the SURMISE output
%
% Returns:
%   - out: Structure storing the simulation's output

data = csvread(filename);
pnum = max(data(:,1)) + 1;
disp([num2str(pnum),' particles found'])
out = cell(1,pnum);
for ip=1:pnum
    out{ip} = struct('xpos',[],'ypos',[],'xvel',[],'yvel',[],'xfrc',[],'yfrc',[]);
end
for il=1:size(data,1)
    out{data(il,1)+1}.xpos(end+1) = data(il,2);
    out{data(il,1)+1}.ypos(end+1) = data(il,3);
    out{data(il,1)+1}.xvel(end+1) = data(il,4);
    out{data(il,1)+1}.yvel(end+1) = data(il,5);
    out{data(il,1)+1}.xfrc(end+1) = data(il,6);
    out{data(il,1)+1}.yfrc(end+1) = data(il,7);
end
end