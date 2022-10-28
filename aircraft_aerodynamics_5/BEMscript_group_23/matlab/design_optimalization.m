function [chord, twist, last_edited] = design_optimalization(last_edited, C_p, previous_C_p, chord, twist, iterationsize, twist_initial)
    if strcmp(last_edited, 'chord') && C_p < previous_C_p
        chord = chord / iterationsize;
        twist = twist + abs(twist_initial/iterationsize);
        last_edited = 'twist';
    elseif strcmp(last_edited, 'twist') && C_p < previous_C_p
        chord = chord * iterationsize;
        twist = twist - abs(twist_initial*iterationsize);
        last_edited = 'chord';
    elseif strcmp(last_edited, 'twist') && C_p > previous_C_p
        twist = twist + abs(twist_initial/iterationsize);
    elseif strcmp(last_edited, 'chord') && C_p > previous_C_p
        chord = chord * iterationsize;
    end
            
