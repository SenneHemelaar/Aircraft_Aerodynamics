%% Function to find the correction for the respective wing tip and root loss.
function Product_of_Froot_Ftip = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, n_blades, axial_induction)
    temp1 = -n_blades/2*(tipradius_R-r_R)/r_R*sqrt( 1+ ((TSR*r_R)^2)/((1-axial_induction)^2));
    Ftip = 2/pi*acos(exp(temp1));
    Ftip(isnan(Ftip)) = 0;
    temp1 = n_blades/2*(rootradius_R-r_R)/r_R*sqrt( 1+ ((TSR*r_R)^2)/((1-axial_induction)^2));
    Froot = 2/pi*acos(exp(temp1));
    Froot(isnan(Froot)) = 0;
    Product_of_Froot_Ftip = Froot*Ftip;