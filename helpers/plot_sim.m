function plot_sim(out)
% PLOT_SIM Plots the evolution of the simulation.

figure
ts = length(out{1}.xpos) %#ok<*NOPRT>
lvl = 3;
xtick = linspace(0,100,1+2^lvl);
ytick = xtick;
for it=1:ts
    axes('XTick',xtick,'YTick',ytick)
    grid on
    hold on
    for ip=1:length(out)
        scatter(out{ip}.xpos(it),out{ip}.ypos(it),200,'.')
        vectarrow([out{ip}.xpos(it),out{ip}.ypos(it)],...
            [out{ip}.xpos(it),out{ip}.ypos(it)]+20.*[out{ip}.xfrc(it),out{ip}.yfrc(it)])
    end
    hold off
    xlim([0 100])
    ylim([0 100])
    axis equal
    title(['it=',num2str(it)])
    drawnow
    pause(0.01)
    clf
end
end