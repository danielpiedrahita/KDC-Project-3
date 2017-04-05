%%% criterion for optimizing initial center of mass location and velocity
function score=criterion(p)
    global markers;
    global markers_world;
    global numFrames;
    global numMarkers; % number of markers

    % pull out parameters
    com = [ p(1) p(2) p(3) ];
    vel = [ p(4) p(5) p(6) ];

    score = 0;

    % calculate distances from com to each marker.
    for j=1:numMarkers
        vector = markers_world(1,(3*(j-1)+1+1):(3*(j-1)+3+1)) - com;
        d1(j) = vector*vector'; %dist^2
    end
    for i=2:numFrames
        for j=1:numMarkers
            vector = markers_world(i,(3*(j-1)+1+1):(3*(j-1)+3+1)) - com - vel*markers_world(i,1);
            dist = vector*vector'; %dist cost
            score = score + (d1(j) - dist)*(d1(j) - dist);
        end
    end
end
