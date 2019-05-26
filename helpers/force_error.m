function [ferrx,ferry]=force_error(out,ip,relative)
% ferr=FORCE_ERROR(out,relative) Compares the simulation's estimation of
% the force acting on a particle compared to the expected (exact) value. It
% takes a cell array produces by READ_OUTPUT as input argument.
%
% Arguments:
%   - out: Cell array containing a SURMISE simulation output.
%   - ip: Index of the considered particle
%   - relative: Optional boolean switch to return a relative force error
%               instead of the absolute one. Default: true.
%
% Returns:
%   - ferr*: Array with the force error of each particle, at each timestep.

narginchk(1,3);
if(nargin==1)
    ip=1;
    relative=true;
elseif(nargin==2)
    relative=true;
end

domain=[0,100]; %Size of the sim's square domain.

% Storage for errors (note: we need to exclude a timestep - the force
% stored at each timestep is the one that WAS applied for the considered
% timestep. Then, this force is related to the positions of the masses at
% the PREVIOUS timestep.)
ferrx = zeros(1,length(out{1}.xpos)-1);
ferry = ferrx;

% For each timestep
for it=2:length(out{1}.xpos)
    % Expected forces on the particle
    fexpx = 0;
    fexpy = 0;
    for ip2=1:length(out)
        if(ip2~=ip && ...
                out{ip}.xpos(it-1) > domain(1) && out{ip}.xpos(it-1) < domain(2) && ...
                out{ip}.ypos(it-1) > domain(1) && out{ip}.ypos(it-1) < domain(2) && ...
                out{ip2}.xpos(it-1) > domain(1) && out{ip2}.xpos(it-1) < domain(2) && ...
                out{ip2}.ypos(it-1) > domain(1) && out{ip2}.ypos(it-1) < domain(2))
            % Note: still miss mass informations, assume one RN
            fexpx = fexpx - out{ip}.mass*out{ip2}.mass/(...
                (out{ip}.xpos(it-1)-out{ip2}.xpos(it-1))^2 + ...
                (out{ip}.ypos(it-1)-out{ip2}.ypos(it-1))^2) * ...
                (out{ip}.xpos(it-1)-out{ip2}.xpos(it-1));
            fexpy = fexpy - out{ip}.mass*out{ip2}.mass/(...
                (out{ip}.xpos(it-1)-out{ip2}.xpos(it-1))^2 + ...
                (out{ip}.ypos(it-1)-out{ip2}.ypos(it-1))^2) * ...
                (out{ip}.ypos(it-1)-out{ip2}.ypos(it-1));
        end
    end
    ferrx(it-1) = out{ip}.xfrc(it)-fexpx;
    ferry(it-1) = out{ip}.yfrc(it)-fexpy;
    if(relative)
        ferrx(it-1) = ferrx(it-1) / fexpx;
        ferry(it-1) = ferry(it-1) / fexpy;
    end
end
end