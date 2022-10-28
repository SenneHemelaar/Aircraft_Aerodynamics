function CT = Glauert_corrected_thrust_coefficient(a)
    CT1 = 1.816;
    CT = CT1 - 4 * (sqrt(CT1)-1) * (1-a);