%ddermadi@stanford.edu
function [traj] = runpathFinder(sessionData, lknn, parameters, s)
    
    data = sessionData;
    traj=[];
    
    try
        
        [G] = pathFinder(data, lknn, parameters, s);
        traj = G.T;
        traj = traj';
       
    catch e
        fprintf('tSPACE failed: %s', getReport(e,'extended'));
        return;        
    end
    
end


