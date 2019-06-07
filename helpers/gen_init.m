% GEN_INIT Generates intial conditions and saves them in an appropriate
% file format that can be readily used to launch a simulation.

%% Basic parameters
domsize = 100;
npart = 100;
dt = 1e-2;
epsilon = 0.1;
theta = 0.4;
maxiter = 1000;
maxwtime = 60;
extratime = 10;
write_freq = 1;
output_path = 'out/';

%% Distribution parameters (custom, user-filled)
sim_idx = 4;
simname = 'default';
savepath= './';
switch sim_idx
    case 0
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1,1);
        simname = 'alluniform';
    case 1
        npart = 100000;
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,1e-7,1e-4);
        simname = 'largeuniform';
    case 2
        npart = 1000000;
        [x,v,m] = ic_alluniform(npart,0,domsize,0,0,0.01,0.1);
        simname = 'vlargeuniform';
    case 3
        npart = 100;
        vrad = 8;
        vstd = 2;
        [x,v,m] = ic_radial(npart,0,domsize,vrad,vstd,0.05,1);
        simname = 'medradial';
    case 4
        npart = 1e6;
        vrad = 1;
        vstd = 0.25;
        clustsz = 15; clustpos=25;
        [x1,v1,m1] = ic_radial(npart/2,0,clustsz,vrad,vstd,0.001,0.1);
        [x2,v2,m2] = ic_radial(npart/2,0,clustsz,vrad,vstd,0.001,0.1);
        x1 = x1 - clustpos;
        x2 = x2 + clustpos;
        x = [x1;x2];
        v = [v1;v2];
        m = [m1;m2];
        clear x1 x2 v1 v2 m1 m2;
        simname = 'clusters';
end
% Expect columnwise storage
data = [x v m];
% m propto volume propto r^2 (2D)
% a=v^2/r and a propto m/r
% v propto sqrt(ra) propto sqrt(m) propto sqrt(r^2)=r

%% Ouputs parameters
fileconf = [savepath,simname,'.conf'];
fileinit = [savepath,simname,'.init'];

%% Generate configuration files
fileID = fopen(fileconf,'w');
fprintf(fileID,['domsize=%.10f\n',...
                'npart=%u\n',...
                'tevol_dt=%.10f\n',...
                'epsilon=%.6f\n',...
                'theta=%.6f\n',...
                'max_iter=%u\n',...
                'walltime=%.10f\n',...
                'extratime=%.10f\n',...
                'write_freq=%u\n',...
                'output_path=%s\n'],...
                domsize,npart,dt,epsilon,theta,maxiter,maxwtime,extratime,...
                write_freq,output_path);
fclose(fileID);
fileID = fopen(fileinit,'w');
fprintf(fileID,'%.25f %.25f %.25f %.25f %.25f\n',data');
fclose(fileID);
