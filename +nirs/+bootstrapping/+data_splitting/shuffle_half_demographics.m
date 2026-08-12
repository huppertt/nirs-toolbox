function data = shuffle_half_demographics(data, relationshipVariable)

    data2 = data;
    tblDemo = nirs.createDemographicsTable(data);

    % Identify the two participant groups
    relationship = string(tblDemo.(relationshipVariable));
    groups = unique(relationship, 'stable');

    if numel(groups) ~= 2
        error('Expected exactly two relationship group');
    end

    lstA = find(relationship == groups(1)); %unchanged group
    lstB = find(relationship == groups(2)); %shuffled group

    nA = length(lstA);
    nB = length(lstB);

    if nA ~= nB
        error('The two relationship groups must have the same number of recordings');
    end

    % Shuffle the group B
    shuffledB = randperm(nB);
    while any(shuffledB == 1:nB)
        shuffledB = randperm(nB);
    end

    % Move the complete demographics between Group B
    % Signal data, time, and probe information remain unchanged
    for i = 1:nB
        data2(lstB(i)).demographics = data(lstB(shuffledB(i))).demographics;
    end

    data = data2;
end