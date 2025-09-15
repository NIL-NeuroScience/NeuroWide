function [Acc] = f_alignAcc(digSona,acc)

trigger = find(diff(digSona) == 1);
sIdx = trigger(1);

delta = round(mean(diff(trigger)));
eIdx = trigger(end)+delta;

Acc = acc(sIdx:eIdx);
%%

% Trigger.num = sum(diff(digSona)>0);
% Trigger.ana = diff(anaSona)>1;
% Trigger.ana(end+1) = 0;
% for i = 1:size(anaSona,1)
%     if Trigger.ana(i) == 1
%         Trigger.ana(i+1:i+3) = 0;
%     end
% end
% 
% Trigger.trigs = find(Trigger.ana);
% diffT = Trigger.trigs(end+1-numchannels)-Trigger.trigs(1);
% diffT = floor(diffT/(Trigger.num/numchannels-1));
% 
% Acc = acc(Trigger.trigs(1):Trigger.trigs(end+1-numchannels)+diffT);

end