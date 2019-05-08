% GEN_INIT Generates intial conditions and saves them in an appropriate
% file format that can be readily used to launch a simulation.

%% Basic parameters
domsize = 100;
npart = 10;
dt = 1e-4;
epsilon = 0.1;
maxiter = 1000;
maxwtime = 60;
extratime = 10;

%% Distribution parameters (custom, user-filled)
sim_idx = 0;
simname = 'default';
savepath= '../input/';
switch sim_idx
    case 0
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'alluniform';
end
% Expect columnwise storage
data = [x v m];

%% Ouputs parameters
fileconf = [savepath,simname,'.conf'];
fileinit = [savepath,simname,'.init'];

%% Generate configuration files
fileID = fopen(fileconf,'w');
fprintf(fileID,['domsize=%.10f\n',...
                'npart=%u\n',...
                'tevol_dt=%.10f\n',...
                'epsilon=%.6f\n',...
                'max_iter=%u\n',...
                'walltime=%.10f\n',...
                'extratime=%.10f\n'],...
                domsize,npart,dt,epsilon,maxiter,maxwtime,extratime);
fclose(fileID);
fileID = fopen(fileinit,'w');
fprintf(fileID,'%.25f %.25f %.25f %.25f %.25f\n',data');
fclose(fileID);