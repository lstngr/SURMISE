% GEN_INIT Generates intial conditions and saves them in an appropriate
% file format that can be readily used to launch a simulation.

%% Basic parameters
domsize = 100;
npart = 10;
dt = 1e-2;
epsilon = 0.1;
maxiter = 1000;
maxwtime = 60;
extratime = 10;

%% Distribution parameters (custom, user-filled)
sim_idx = 3;
simname = 'default';
savepath= './';
switch sim_idx
    case 0
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'alluniform';
    case 1
        npart = 100000;
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'largeuniform';
    case 2
        npart = 100;
        vrad = 8;
        vstd = 2;
        [x,v,m] = ic_radial(npart,0,domsize,vrad,vstd,0.05,1);
        simname = 'medradial';
    case 3
        npart=1;
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'easy';
    case 4
        npart=2;
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'easytoo';
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