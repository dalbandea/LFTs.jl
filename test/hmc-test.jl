

function reversibility!(lftws::LFTs.AbstractLFT, hmcws::LFTs.AbstractHMC)
    LFTs.molecular_dynamics!(lftws, hmcws)
    LFTs.flip_momenta_sign!(hmcws)
    LFTs.molecular_dynamics!(lftws, hmcws)
    return nothing
end

function force_test(lftws::LFTs.AbstractLFT, hmcws::LFTs.AbstractHMC, epsilon)
    lftws2 = deepcopy(lftws)

    Si = action(lftws, hmcws)

    F_ana = analytic_force(lftws, hmcws)

    fld = get_field(lftws2)

    F_diff = 0.0

    for i in 1:length(fld)
        fld[i] = infinitesimal_transformation(fld[i], epsilon, lftws2)

        # Final action
        Sf = action(lftws2, hmcws)

        # Numerical force at point i
        F_num = (Sf - Si)/epsilon

        # Difference
        F_diff += F_ana[i] + F_num

        fld[i] = infinitesimal_transformation(fld[i], -epsilon, lftws2)
    end

    return F_diff / length(fld)
end
