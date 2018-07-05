function myassignin(str,V)
    assignin('base','assignin_temp',V);
    evalin('base',[str ' = assignin_temp;']);
    evalin('base','clear assignin_temp');