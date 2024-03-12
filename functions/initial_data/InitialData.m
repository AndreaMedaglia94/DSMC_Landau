function [vx, vy, vz] = InitialData(p_sch, p_sim)

fprintf('Initializing the particles: ');

if strcmp(p_sch.test, 'BKW')
    [vx, vy, vz] = InitialDataBKW(p_sim);
elseif strcmp(p_sch.test, 'Trub')
    [vx, vy, vz] = InitialDataTrubnikov(p_sim);
end

fprintf('Done \n\n');

end