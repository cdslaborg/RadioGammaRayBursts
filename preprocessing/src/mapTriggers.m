function Map = mapTriggers(TrigTable1,TrigRadio,FreqRadio)
    lenTrigTable1 = length(TrigTable1);
    lenTrigRadio = length(TrigRadio);
    maxLen = max(lenTrigTable1,lenTrigRadio);
    Map.Trig = cell(maxLen,1);
    Map.IndxTable1 = zeros(maxLen,1);
    Map.IndxRadio = zeros(maxLen,1);
    counter = 0;
    for i = 1:lenTrigTable1
        j = 0;
        while j<lenTrigRadio
            j = j + 1;
            if strcmp(TrigTable1{i},TrigRadio{j}) && FreqRadio(j)==8.46
                disp([TrigTable1(i),TrigRadio(j)])
                counter = counter + 1;
                Map.Trig{counter} = TrigTable1{i};
                Map.IndxTable1(counter) = i;
                Map.IndxRadio(counter) = j;
                %break;
                j = lenTrigRadio + 1;
            end
        end
    end
    Map.Trig = Map.Trig(1:counter);
    Map.IndxTable1 = Map.IndxTable1(1:counter);
    Map.IndxRadio = Map.IndxRadio(1:counter);
    Map.count = counter;
end