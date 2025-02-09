using QuadGK
using NLsolve
include("constants.jl")
println("Nf_o, Nc_o, Lambda_o, G, m_o are global variables.")

function np(p, T, mu, M)
    """
    Computes the Fermi-Dirac occupation number for quarks.

    Args:
        p: The momentum of the quark.
        T: Temperature (in appropriate energy units, e.g., MeV).
        mu: Chemical potential (MeV).
        M: Effective mass of the quark (MeV).
    
    Returns:
        A Float64 representing the quark occupation probability.
    """
    Ep = sqrt(p^2 + M^2)
    return 1.0/(exp((Ep-mu)/T) + 1.0)
end

function np_bar(p, T, mu, M)
    """
    Computes the Fermi-Dirac occupation number for antiquarks.

    Args:
        p: The momentum of the antiquark.
        T: Temperature (MeV).
        mu: Chemical potential (MeV).
        M: Effective mass of the quark (MeV).
    
    Returns:
        A Float64 representing the antiquark occupation probability.
    """
    Ep = sqrt(p^2 + M^2)
    return 1.0/(exp((Ep+mu)/T) + 1.0)
end

function num_dens(T, mu, M)
    """
    Computes the total quark number density.

    Args:
        T: Temperature (MeV).
        mu: Chemical potential (MeV).
        M: Effective quark mass (MeV).
    
    Returns:
        The quark number density (with appropriate prefactors), as a Float64.
    """
    function integrand_numdens(p)
        if T > 0.0
            return (np(p, T, mu, M) - np_bar(p, T, mu, M))*p^2  # polar coordinates!
        else
            Ep = sqrt(p^2 + M^2)
            step_function = x -> x >= 0 ? 1.0 : 0.0
            return step_function(mu-Ep)*p^2
        end
    end
    integral, _ = quadgk(integrand_numdens, 0.0, Lambda_o)
    return 4.0 * Nf_o * Nc_o / (2.0 * π)^2 * integral  # 4π is multiplied ! 
end

function rhoB_o(T,μ)
    """
    Computes the normalized baryon number density.
    
    Args:
        T: Temperature (MeV).
        μ: Chemical potential (MeV).
    
    Returns:
        The normalized baryon density.
    """
    M_sol = find_M(T, μ)
    return num_dens(T, μ, M_sol)/(hbarc^3)/3.0/rho0
end


function gap_eq(M; T=0.0, mu=0.0)
    """
    Computes the gap equation (the self-consistency condition) for the quark mass.
    
    The gap equation is defined via an integral over momentum space. For T > 0, the integrand uses the 
    Fermi-Dirac distributions for quarks and antiquarks. For T = 0, a step function is used instead.
    The gap equation is given by:
    
        gap_eq(M) = 8.0 * Nf_o * Nc_o * G/(2π)^2 * ∫ dp [integrand] + m_o - M = 0.0
    
    Args:
        M: The effective quark mass (MeV) to be solved for.
        T (keyword): Temperature (MeV), default 0.
        mu (keyword): Chemical potential (MeV), default 0.
    
    Returns:
        The value of the gap equation. The correct M should satisfy gap_eq(M) = 0.
    """
    function integrand_gapeq(p)
        if T > 0.0
            Ep = sqrt(p^2 + M^2)
            return M/Ep*(1.0 - np(p, T, mu, M) - np_bar(p, T, mu, M))*p^2  # polar coordinates!
        else
            Ep = sqrt(p^2 + M^2)
            step_function = x -> x >= 0 ? 1.0 : 0.0
            return M/Ep*(1.0 - step_function(mu-Ep))*p^2
        end
    end
    integral, _ = quadgk(integrand_gapeq, 0.0, Lambda_o)
    return 8.0*Nf_o*Nc_o*G/(2.0*π)^2*integral + m_o - M # 4π is multiplied ! 
end

function find_M(Temp, μ)
    """
    Solves the gap equation numerically to find the effective quark mass M.
    
    This function uses a simple iterative approach combined with the nonlinear solver `nlsolve` to
    find the value of M that satisfies gap_eq(M; T, mu) = 0 for a given temperature (Temp) and chemical
    potential (μ).
    
    Args:
        Temp: Temperature (MeV).
        μ: Chemical potential (MeV).
    
    Returns:
        The effective quark mass M (MeV) that satisfies the gap equation.
    """
    M_guess = 1e3
    M_temp = M_guess
    m_step = 1e2
    while abs(gap_eq(M_temp; T=Temp, mu=μ)) > 1e-4 && M_guess > -1e-4
        M_guess -= m_step
        sol = nlsolve(x -> [gap_eq(x[1]; T=Temp, mu=μ)], [M_guess])
        M_temp = sol.zero[1]
    end
    return M_temp  # Return the solution for M
end