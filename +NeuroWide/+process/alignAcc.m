function Acc = alignAcc(digSona,acc)

trigger = find(diff(digSona) == 1);
sIdx = trigger(1);

delta = round(mean(diff(trigger)));
eIdx = trigger(end)+delta;

Acc = acc(sIdx:eIdx);

end