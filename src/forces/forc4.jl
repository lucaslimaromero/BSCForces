using SpecialFunctions
using Plots
using LegendrePolynomials

# Struct de parâmetros físicos
struct ParametrosSimulacao
    f::Float64; c₀::Float64; ρ₀::Float64; P₀::Float64
    a::Float64; ρ_s::Float64; cL_s::Float64; cT_s::Float64
    ν::Int; α::Float64
    λ::Float64; k::Float64; ω::Float64; ϕ₀_amp::Float64
end

function ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_deg, a_div_λ)
    λ = c₀ / f
    a = a_div_λ * λ
    k = 2π / λ
    ω = 2π * f
    ϕ₀_amp = P₀ / (ω * ρ₀)
    α = deg2rad(α_deg)
    return ParametrosSimulacao(f, c₀, ρ₀, P₀, a, ρ_s, cL_s, cT_s, ν, α, λ, k, ω, ϕ₀_amp)
end

# BSCs deslocados (completos)
function bsccalc(n, m, p::ParametrosSimulacao, ρ₀, ϕ₀)
    if n < abs(m); return 0.0 + 0.0im; end
    kρ = p.k * sin(p.α)
    logfac = lgamma(n - abs(m) + 1) - lgamma(n + abs(m) + 1)
    Plmval = Plm(cos(p.α), n, abs(m))
    bessel = besselj(m - p.ν, kρ * ρ₀)
    fase = cis(-(m - p.ν) * ϕ₀)
    return im^(n - m) * (2n + 1) * exp(logfac) * Plmval * bessel * fase
end

# Coeficientes de espalhamento Rn
function calcular_Rn(n, p::ParametrosSimulacao)
    ω, k₀, a = p.ω, p.k, p.a
    kL = ω / p.cL_s
    kT = ω / p.cT_s
    x, xL, xT = k₀ * a, kL * a, kT * a

    ξ(n, z, f) = Float64(n) - z * f(n + 1, z) / f(n, z)
    sphericalhankelh(n, z) = sphericalbesselj(n, z) + im * sphericalbessely(n, z)

    ξj = ξ(n, x, sphericalbesselj)
    ξh = ξ(n, x, sphericalhankelh)
    ξL = ξ(n, xL, sphericalbesselj)
    ξT = ξ(n, xT, sphericalbesselj)

    dens_ratio = p.ρ_s / p.ρ₀

    num1 = 2n * (n + 1)
    den1 = n^2 + n - 0.5xT^2 + 2 * (n^2 - 1) / ((ξT - 1) / (ξL - 1) - 1)
    num = ξj * (dens_ratio * ξL + num1 / den1) - n*(n+1)*dens_ratio
    den = ξh * (dens_ratio * ξL + num1 / den1) - n*(n+1)*dens_ratio

    return abs(den) < 1e-30 ? 0.0 + 0.0im : -num / den
end

# Força axial Eq. (6c)
function calcular_forca_axial(p::ParametrosSimulacao, n_max, bscs_2D, α_coeffs, β_coeffs)
    soma = 0.0
    for n in 0:n_max-1
        αn, βn = α_coeffs[n+1], β_coeffs[n+1]
        αn1, βn1 = α_coeffs[n+2], β_coeffs[n+2]
        D1 = αn + αn1 + 2*(αn*αn1 + βn*βn1)
        D2 = βn1 - βn + 2*(βn1*αn - αn1*βn)
        for m in -n:n
            Qnm = 2 * exp(lgamma(n + m + 1) - lgamma(n - m + 1) - log(2n + 1) - log(2n + 3))
            Anm = bscs_2D[n+1, m + n_max + 1]
            An1m = bscs_2D[n+2, m + n_max + 1]
            soma += 2*(n + m + 1)*Qnm * (real(Anm * conj(An1m)) * D2 - imag(Anm * conj(An1m)) * D1)
        end
    end
    return -(π * p.ρ₀ * p.ϕ₀_amp^2 / 2) * soma
end

# Plot geral
function gerar_figura4()
    f = 1e6; P₀ = 0.1e6; c₀ = 1500.0; ρ₀ = 1000.0
    cL_s = 2350.0; cT_s = 1120.0; ρ_s = 1050.0
    ν = 1; n_max = 50
    a_div_λ_range = 0.001:0.005:0.5
    α_vals = [50.0, 60.0, 70.0]

    plt = plot(xlabel="a/λ", ylabel="F_z [N]",
        title="Força Axial em Esfera de Poliestireno (ν=1)", legend=:topleft)

    for α_deg in α_vals
        Fz_array = zeros(length(a_div_λ_range))
        for (i, a_div_λ) in enumerate(a_div_λ_range)
            p = ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_deg, a_div_λ)
            α_coeffs = zeros(n_max + 1); β_coeffs = zeros(n_max + 1)
            for n in 0:n_max
                Rn = calcular_Rn(n, p)
                α_coeffs[n+1] = real(Rn)
                β_coeffs[n+1] = imag(Rn)
            end
            bscs_2D = zeros(ComplexF64, n_max + 2, 2n_max + 1)
            for n in 0:n_max+1
                for m in -n:n
                    col = m + n_max + 1
                    if col < 1 || col > size(bscs_2D, 2)
                        continue
                    end
                    bscs_2D[n+1, col] = bsccalc(n, m, p, 0.0, 0.0)
                end
            end
            Fz_array[i] = calcular_forca_axial(p, n_max, bscs_2D, α_coeffs, β_coeffs)
        end
        plot!(plt, a_div_λ_range, Fz_array, label="β = $(α_deg)°", lw=2)
    end

    savefig(plt, "figura4_correta.png")
    println("Figura 4 salva como 'figura4_correta.png'")
end

# Execute
gerar_figura4()
