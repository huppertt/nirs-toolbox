function data2 = randomize_hyperscanning_pairs(data, linkVariable)

    tblDemo = nirs.createDemographicsTable(data);

    lstlinkVariable = unique(tblDemo(:,linkVariable),'rows');

    lstA = find(ismember(tblDemo.relationship, lstlinkVariable.relationship(1)));
    dataA = data(lstA);
    
    lstB = find(ismember(tblDemo.relationship, lstlinkVariable.relationship(2)));
    dataB = data(lstB);

    nA = length(dataA);
    nB = length(dataB);

    if nA ~= nB
        error('dataA and dataB must contain the same number of recordings');
    end

    % Shuffle B recordings, ensuring no original A-B pair is retained
    shuffledB = randperm(nB);
    while any(shuffledB == 1:nA)
        shuffledB = randperm(nB);
    end

    for i = 1:nA
        recordingA = dataA(i);
        recordingB = dataB(shuffledB(i));

        % Trim each pair to their shared available time length
        nTime = min(length(recordingA.time), length(recordingB.time));

        pairedData(i).subjectA = recordingA;
        pairedData(i).subjectA.time = recordingA.time(1:nTime);
        pairedData(i).subjectA.data = recordingA.data(1:nTime,:);

        pairedData(i).subjectB = recordingB;
        pairedData(i).subjectB.time = recordingB.time(1:nTime);
        pairedData(i).subjectB.data = recordingB.data(1:nTime,:);

        % Keep the record of the shuffled pairing
        pairedData(i).subjectA_index = i;
        pairedData(i).subjectB_index = shuffledB(i);
    end
    
    data2 = reshape([pairedData.subjectA; pairedData.subjectB], [], 1);

end