using QuadGK  # 数値積分のためのパッケージ
using NLsolve

include("constants.jl")
println("Nf, Nc, Lambda, Gs, Gv, m are global variables.")

function np(p, T, mu, M)
    """
    Computes the Fermi occupation number for quarks.

    Args:
        p: Momentum.
        T: Temperature.
        mu: Chemical potential.
        M: Effective mass.

    Returns:
        Occupation probability for quarks.
    """
    Ep = sqrt(p^2 + M^2)
    return 1.0/(exp((Ep-mu)/T) + 1.0)
end

function np_bar(p, T, mu, M)
    """
    Computes the Fermi occupation number for antiquarks.

    Args:
        p: Momentum.
        T: Temperature.
        mu: Chemical potential.
        M: Effective mass.

    Returns:
        Occupation probability for antiquarks.
    """
    Ep = sqrt(p^2 + M^2)
    return 1.0/(exp((Ep+mu)/T) + 1.0)
end

function OmegaM(Temp, μ_t, Mass)
    """
    Computes the free Fermi-gasthermodynamic potential (grand potential) for given 
    temperature, chemical potential, and mass.

    Args:
        Temp: Temperature.
        μ_t: Effective chemical potential.
        Mass: Effective mass.

    Returns:
        The thermodynamic potential OmegaM.
    """
    function integrand_OmegaM(p)
        Ep = sqrt(p^2+Mass^2)
        if Temp > 0.0
            return p^2 * (Ep + Temp * log(1.0 + exp(-(Ep - μ_t)/Temp)) + Temp*log(1.0 + exp(-(Ep + μ_t)/Temp)))
        else
            step_function = x -> x >= 0 ? 1.0 : 0.0
            return p^2 * (Ep + (μ_t-Ep)*step_function(μ_t-Ep))
        end
    end
    integral, _ = quadgk(integrand_OmegaM, 0.0, Lambda)  # Regularized integral
    return -4.0 * Nf * Nc / (2.0*π)^2 * integral

end

function Omega_temp(T, μ, M, μ_tilde)
    """
    Computes the total thermodynamic potential including interaction terms.

    Args:
        T: Temperature.
        μ: Chemical potential.
        M: Effective mass.
        μ_tilde: Effective chemical potential.

    Returns:
        The total Omega value.
    """
    return OmegaM(T, μ_tilde, M) + (M-m)^2/(4.0*Gs) - (μ - μ_tilde)^2/(4.0*Gv) +1.9397179193981308e10 # Any const can be added
end

function dOmegadM(T, μ, Mass, μ_tilde)
    """
    Computes the derivative of Omega with respect to mass M.

    Args:
        T: Temperature.
        μ: Chemical potential.
        Mass: Effective mass.
        μ_tilde: Effective chemical potential.

    Returns:
        The derivative value.
    """
    function integrand_dOegadM(p)
        Ep = sqrt(p^2+Mass^2)
        if T > 0.0
            return p^2 * Mass/Ep*(1-np(p,T,μ_tilde, Mass)-np_bar(p,T,μ_tilde, Mass))
        else
            step_function = x -> x >= 0 ? 1.0 : 0.0
            return p^2 * Mass/Ep*(1-step_function(μ_tilde-Ep))
        end
            
    end
    integral, _ = quadgk(integrand_dOegadM, 0.0, Lambda)
    return (Mass-m)/(2.0*Gs) -  4.0*Nf*Nc/(2.0*π)^2 *integral
end

function dOmegadmutilde(T, μ, Mass, μ_tilde)
    """
    Computes the derivative of Omega with respect to effective chemical potential μ_tilde.

    Args:
        T: Temperature.
        μ: Chemical potential.
        Mass: Effective mass.
        μ_tilde: Effective chemical potential.

    Returns:
        The derivative value.
    """
    function integrand_dOmegadmutilde(p)
        if T > 0.0
            return  p^2 * (np(p, T, μ_tilde, Mass)-np_bar(p, T, μ_tilde, Mass))
        else
            Ep = sqrt(p^2+Mass^2)
            step_function = x -> x >= 0 ? 1.0 : 0.0
            return p^2 * step_function(μ_tilde - Ep)
        end
    end
    integral, _ = quadgk(integrand_dOmegadmutilde, 0.0, Lambda)
    return (μ - μ_tilde)/(2.0*Gv) - 4.0*Nf*Nc/(2.0*π)^2 * integral
end

function solve_system(T, μ)
    """
    Solves the system of equations dOmega/dM = 0 and dOmega/dμ_tilde = 0.

    Args:
        T: Temperature.
        μ: Chemical potential.

    Returns:
        The values of (M, μ_tilde) that minimize Omega.
    """
    M_guess = [4e2, m]
    μ_guess = [μ, μ/2.0]

    f!(F, x) = begin
        F[1] = dOmegadM(T, μ, x[1], x[2])  # First equation
        F[2] = dOmegadmutilde(T, μ, x[1], x[2])   # Second equation
    end

    sol1 =  nlsolve(f!, [M_guess[1], μ_guess[1]])
    sol2 =  nlsolve(f!, [M_guess[2], μ_guess[2]])
    sol1_temp = sol1.zero
    sol2_temp = sol2.zero

    if Omega_temp(T, μ, sol1_temp[1], sol1_temp[2]) <= Omega_temp(T, μ, sol2_temp[1], sol2_temp[2])
        # println(sol1_temp)
        return sol1_temp  # Returns [M, μ_tilde]
    else
        return sol2_temp
    end
end



function Omega(T, μ)
    """
    Computes the grand potential Omega(T, μ).

    Args:
        T: Temperature.
        μ: Chemical potential.

    Returns:
        Omega relative to its vacuum value.
    """
    M_sol, μ_tilde_sol = solve_system(T, μ)
    Mvac, mu_tilde_vac = solve_system(0.0, 0.0)
    return Omega_temp(T, μ, M_sol, μ_tilde_sol) - Omega_temp(0.0, 0.0, Mvac, 0.0)
    # To guarantee p > 0.0, eps(Float64 is added.
end

function num_dens(T, mu, M)
    """
    Computes the quark number density.

    Args:
        T: Temperature.
        mu: Chemical potential.
        M: Effective mass.

    Returns:
        The quark number density.
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
    integral, _ = quadgk(integrand_numdens, 0.0, Lambda)
    return 4.0 * Nf * Nc / (2.0 * π)^2 * integral  # 4π is multiplied ! 
end

function rhoB(T,μ)
    """
    Computes the normalized baryon number density for given temperature T and chemical potential μ.
      
    Args:
        T::Float64: Temperature (MeV).
        μ::Float64: Quark chemical potential (MeV).
        
    Returns:
        The normalized baryon number density (dimensionless or in units of saturation density, depending on rho0).
    """
    M_sol, μ_tilde_sol = solve_system(T, μ)
    return num_dens(T, μ_tilde_sol, M_sol)/(hbarc^3)/3.0/rho0
end

function pressure(T, μ)
    """
    Computes the pressure.

    Args:
        T: Temperature.
        μ: Chemical potential.

    Returns:
        Pressure value.
    """
    M_sol, μ_tilde_sol = solve_system(T, μ)
    Mvac, mu_tilde_vac = solve_system(0.0, 0.0)
    return - (Omega_temp(T, μ, M_sol, μ_tilde_sol) - Omega_temp(0.0, 0.0, Mvac, 0.0))
    # To guarantee p > 0.0, eps(Float64 is added.
end

function energy_dens(T, μ)
    """
    Computes the energy density.

    Args:
        T: Temperature.
        μ: Chemical potential.

    Returns:
        Energy density value.
    """
    M_sol, μ_tilde_sol = solve_system(T, μ)
    if T > 0.0
        h = 1e-8
        s = -(Omega_temp(T+h, μ, M_sol, μ_tilde_sol) - Omega_temp(T, μ, M_sol, μ_tilde_sol))/h
    else
        s = 0.0
    end
    return -pressure(T, μ) + T*s + μ*num_dens(T, μ_tilde_sol, M_sol)
end