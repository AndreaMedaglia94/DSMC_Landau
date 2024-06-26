function p_sch = set_scheme_parameters

% % % choose the test
p_sch.test = 'BKW' ;
% p_sch.test = 'Trub' ;

% % % choose the DSMC scheme: Nanbu-Babovski (NB)
p_sch.coll = 'NB';

% % % choose the time order
p_sch.timeorder = 2 ;

% % % choose the approximated surrogate kernel D^i_* i=2,3
% p_sch.kernel = 'D1' ;
% p_sch.kernel = 'D2' ;
p_sch.kernel = 'D3' ;

% % % choose the collision scenario
% p_sch.pot = 'Coulomb' ;
p_sch.pot = 'Maxwell' ;

% % % plot
p_sch.init_conditions = 'NO' ;
p_sch.transient       = 'NO' ;

if strcmp(p_sch.test, 'BKW')
    if strcmp(p_sch.pot, 'Coulomb')
        fprintf('Exact BKW solution works only with Maxwellian particles! \n');
        stop
    end
    fprintf('BKW Test \n');
elseif strcmp(p_sch.test, 'Trub')
    fprintf('Trubnikov Test \n');
else
    fprintf('Error: no test has been selected \n');
    stop    
end


if strcmp(p_sch.coll, 'NB')
    fprintf('Nanbu-Babovski scheme \n');
else
    fprintf('Error: no scheme for the collisional operator has been selected \n');
    stop     
end

if strcmp(p_sch.kernel, 'D1')
    fprintf('Approximated collisional kernel D1 \n');
elseif strcmp(p_sch.kernel, 'D2')
    fprintf('Approximated collisional kernel D2 \n');
elseif strcmp(p_sch.kernel, 'D3')
    fprintf('Approximated collisional kernel D3 \n');    
else
    fprintf('Error: no approximated collisional kernel has been selected \n');
    stop     
end

if strcmp(p_sch.pot, 'Maxwell')
    fprintf('Interaction potential: Maxwell \n');
elseif strcmp(p_sch.pot, 'Coulomb')
    fprintf('Interaction potential: Coulomb \n');
else
    fprintf('Error: no interaction potential has been selected \n');
    stop     
end

if p_sch.timeorder == 1
    fprintf('First order Forward Eueler scheme in time \n');
elseif p_sch.timeorder == 2
    fprintf('Second order Midpoint method in time \n');
else
    fprintf('Error: no time scheme has been selected \n');
    stop     
end

end