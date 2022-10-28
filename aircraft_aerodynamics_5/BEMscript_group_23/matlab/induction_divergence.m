% Step 7 of BEM iterations. Inputs (a, a_last, a_prime, a_prime_last,
% error_limit). Returns boolean true if converged

function res = induction_divergence(a, a_last, a_prime, a_prime_last, error_limit)

%Outputs a boolean true or false (1/0). Error limit can be reduced to
%increase accuracy of BEM

a_error = abs(a - a_last); %using absolute value, finding errors to compare with error_limit
a_prime_error = abs(a_prime - a_prime_last);

if (a_error >= error_limit) && (a_prime_error >= error_limit)
    res = true; %Yes, diverged
else 
    res = false; %No, did not diverge (it converged. Good)
    
end %if

end %fun