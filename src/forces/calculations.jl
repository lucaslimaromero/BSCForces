# calculations.jl

using SpecialFunctions
using LegendrePolynomials

# ==============================================================================
# ESTRUTURA PARA ORGANIZAR OS PARÂMETROS (Mantida)
# ==============================================================================
struct ParametrosSimulacao
    f::Float64; c₀::Float64; ρ₀::Float64; P₀::Float64
    a::Float64; ρ_s::Float64; cL_s::Float64; cT_s::Float64
    ν::Int; α_axicon::Float64
    λ::Float64; k::Float64; ω::Float64; ϕ₀_amp::Float64
end

function ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_axicon_deg, a_div_λ)
    λ = c₀ / f
    a = a_div_λ * λ
    k = 2π / λ
    ω = 2π * f
    ϕ₀_amp = P₀ / (ω * ρ₀)
    α_axicon_rad = deg2rad(α_axicon_deg)
    return ParametrosSimulacao(f, c₀, ρ₀, P₀, a, ρ_s, cL_s, cT_s, ν, α_axicon_rad, λ, k, ω, ϕ₀_amp)
end

# calculations.jl

using SpecialFunctions
using LegendrePolynomials

# ... (Mantenha a struct e as funções de derivada como estão) ...
# ... (Mantenha a função 'calcular_coeficientes_espalhamento' que funcionou para a Fig. 4) ...
# ... (Mantenha as funções 'bsccalc' e 'calcular_forcas') ...


# ==============================================================================
# CÁLCULO DOS COEFICIENTES DE ESPALHAMENTO (APROXIMAÇÃO DE RAYLEIGH)
# Para a << λ, esta é a forma mais simples e correta de calcular os coeficientes.
# Esta função será usada especificamente para a Figura 2.
# ==============================================================================
# ==============================================================================
# CÁLCULO DOS COEFICIENTES DE ESPALHAMENTO (APROXIMAÇÃO DE RAYLEIGH)
# CORREÇÃO FINAL DE SINAL
# ==============================================================================
function calcular_coeffs_rayleigh(n_max::Int, p::ParametrosSimulacao)
    α_coeffs = zeros(Float64, n_max + 2)
    β_coeffs = zeros(Float64, n_max + 2)
    x = p.k * p.a
    
    K_fluido = p.ρ₀ * p.c₀^2
    K_solido = p.ρ_s * p.cL_s^2
    
    f₀ = 1.0 - (K_fluido / K_solido)
    f₁ = (2.0 * (p.ρ_s / p.ρ₀ - 1.0)) / (2.0 * p.ρ_s / p.ρ₀ + 1.0)
    
    # CORREÇÃO DE SINAL: Trocamos -im por +im para alinhar com a convenção da força.
    if n_max >= 0
        R0 = im * (x^3 / 3.0) * f₀
        α_coeffs[1] = real(R0)
        β_coeffs[1] = imag(R0)
    end
    if n_max >= 1
        R1 = im * (x^3 / 3.0) * f₁
        α_coeffs[2] = real(R1)
        β_coeffs[2] = imag(R1)
    end
    
    return α_coeffs, β_coeffs
end

# ==============================================================================
# FUNÇÕES DE DERIVADA (Mantidas)
# ==============================================================================
function sphericalbesselj_prime(n, z)
    return n == 0 ? -sphericalbesselj(1, z) : sphericalbesselj(n - 1, z) - ((n + 1) / z) * sphericalbesselj(n, z)
end
function sphericalbessely_prime(n, z)
    return n == 0 ? -sphericalbessely(1, z) : sphericalbessely(n - 1, z) - ((n + 1) / z) * sphericalbessely(n, z)
end
function sphericalhankelh_prime(n, z)
    return sphericalbesselj_prime(n, z) + im * sphericalbessely_prime(n, z)
end

# ==============================================================================
# CÁLCULO DOS COEFICIENTES DE ESPALHAMENTO (VERSÃO FINAL - HICKLING/FARAN)
# ==============================================================================
function calcular_coeficientes_espalhamento(n_max::Int, p::ParametrosSimulacao)
    x = p.k * p.a
    xL = (p.ω / p.cL_s) * p.a
    xT = (p.ω / p.cT_s) * p.a
    ρ_ratio = p.ρ_s / p.ρ₀
    c_ratio_L = p.cL_s / p.c₀
    c_ratio_T = p.cT_s / p.c₀

    α_coeffs = zeros(Float64, n_max + 2)
    β_coeffs = zeros(Float64, n_max + 2)

    for n in 0:n_max + 1
        # Funções de Bessel para o fluido
        jn_x = sphericalbesselj(n, x)
        yn_x = sphericalbessely(n, x)
        jnp_x = sphericalbesselj_prime(n, x)
        ynp_x = sphericalbessely_prime(n, x)

        # Log-derivadas para a esfera
        # D_n(z) = z * j_n'(z) / j_n(z)
        Dn_xL = xL * sphericalbesselj_prime(n, xL) / sphericalbesselj(n, xL)
        Dn_xT = xT * sphericalbesselj_prime(n, xT) / sphericalbesselj(n, xT)

        # Termo de impedância da esfera Zₙ (extremamente sensível)
        term1 = (2 * n * (n + 1) * Dn_xL) / ((n^2 + n - 0.5 * xT^2) * Dn_xL - (n^2 - 1) * n * (n + 1))
        term2 = (Dn_xT) / (Dn_xT - 1)
        Z_n_shear = (n^2 + n - 0.5 * xT^2 - n * (n + 1) * term2) * term1
        
        # Impedância total da esfera
        Z_n_total = -ρ_ratio * (c_ratio_L^2 * Dn_xL + c_ratio_T^2 * Z_n_shear)
        
        # Coeficiente de espalhamento Rₙ (ou Cₙ em alguns textos)
        numerador_Rn = (x * jnp_x / jn_x) * jn_x + Z_n_total * jn_x
        denominador_Rn = (x * (jnp_x + im*ynp_x) / (jn_x + im*yn_x)) * (jn_x + im*yn_x) + Z_n_total * (jn_x + im*yn_x)
        
        Rn = -numerador_Rn / denominador_Rn

        if isnan(Rn); Rn=0.0+0.0im; end

        α_coeffs[n+1] = real(Rn)
        β_coeffs[n+1] = imag(Rn)
    end
    return α_coeffs, β_coeffs
end


# ==============================================================================
# BSC E CÁLCULO DE FORÇAS (Mantidos como antes)
# ==============================================================================
# function bsccalc(n::Int, m::Int, p::ParametrosSimulacao; ρ_desloc=0.0, ϕ_desloc=0.0)
#     if abs(m) > n; return ComplexF64(0.0); end
#     m_abs = abs(m)
#     Pnm_val = Plm(cos(p.α_axicon), n, m_abs)
#     log_fact_ratio = lgamma(n - m_abs + 1) - lgamma(n + m_abs + 1)
#     A_on_axis = im^(n - m) * (2n + 1) * Pnm_val
#     if m < 0
#         A_on_axis *= (-1)^m_abs * exp(-log_fact_ratio)
#     else
#         A_on_axis *= exp(log_fact_ratio)
#     end
#     if ρ_desloc == 0.0
#         return (m == p.ν) ? A_on_axis : ComplexF64(0.0)
#     else
#         k_rho = p.k * sin(p.α_axicon)
#         termo_bessel = besselj(m - p.ν, k_rho * ρ_desloc)
#         termo_fase = cis(-(m - p.ν) * ϕ_desloc)
#         return A_on_axis * termo_bessel * termo_fase
#     end
# end
function bsccalc(n::Int, m::Int, p::ParametrosSimulacao; ρ_desloc=0.0, ϕ_desloc=0.0)
    # Condição fundamental: |m| <= n
    if abs(m) > n
        return ComplexF64(0.0)
    end
    
    # --- Cálculo do coeficiente no eixo (On-axis BSC) ---
    m_abs = abs(m)
    Pnm_val = Plm(cos(p.α_axicon), n, m_abs)
    
    A_on_axis = 0.0 + 0.0im
    
    # A fórmula para m<0 é diferente da de m>=0 devido às convenções
    # dos polinômios de Legendre. Esta é a implementação correta.
    if m >= 0
        log_fact_ratio = lgamma(n - m + 1) - lgamma(n + m + 1)
        A_on_axis = im^(n - m) * (2n + 1) * exp(log_fact_ratio) * Pnm_val
    else # m < 0
        # A fórmula derivada da relação de simetria para m negativo
        A_on_axis = (-1)^m_abs * im^(n + m_abs) * (2n + 1) * Pnm_val
    end

    # --- Aplicação do deslocamento (Off-axis) ---
    # Esta parte do código já estava correta.
    if ρ_desloc == 0.0
        # No eixo, só m=ν contribui para o feixe de Bessel puro
        return (m == p.ν) ? A_on_axis : ComplexF64(0.0)
    else
        k_rho = p.k * sin(p.α_axicon)
        termo_bessel = besselj(m - p.ν, k_rho * ρ_desloc)
        termo_fase = cis(-(m - p.ν) * ϕ_desloc)
        return A_on_axis * termo_bessel * termo_fase
    end
end

function calcular_forcas(p::ParametrosSimulacao, n_max::Int, bscs_2D, α_coeffs, β_coeffs)
    soma_Fx, soma_Fy, soma_Fz = 0.0, 0.0, 0.0
    for n in 0:n_max
        αn, βn = α_coeffs[n+1], β_coeffs[n+1]
        αn_p1, βn_p1 = α_coeffs[n+2], β_coeffs[n+2]
        D1 = αn + αn_p1 + 2 * (αn * αn_p1 + βn * βn_p1)
        D2 = βn_p1 - βn + 2 * (βn_p1 * αn - αn_p1 * βn)
        for m in -n:n
            log_Q_nm = log(2) + lgamma(n + m + 1) - (log(2n + 1) + log(2n + 3) + lgamma(n - m + 1))
            Q_nm = exp(log_Q_nm)
            Anm = bscs_2D[n + 1, m + n_max + 2]
            if n < n_max
                Anp1_m = bscs_2D[n + 2, m + n_max + 2]
                termo_acopl_z = Anm * conj(Anp1_m)
                soma_Fz += 2 * (n + m + 1) * Q_nm * (real(termo_acopl_z) * D2 - imag(termo_acopl_z) * D1)
            end
            V_nm = Float64((n + m + 1) * (n + m + 2))
            Anp1_mp1 = (abs(m + 1) <= n + 1) ? bscs_2D[n + 2, (m + 1) + n_max + 2] : 0.0
            Anp1_mm1 = (abs(m - 1) <= n + 1) ? bscs_2D[n + 2, (m - 1) + n_max + 2] : 0.0
            termo_acopl_p1 = Anm * conj(Anp1_mp1)
            termo_acopl_m1 = Anm * conj(Anp1_mm1)
            Im_term_fx = V_nm * imag(termo_acopl_p1) - imag(termo_acopl_m1)
            Re_term_fx = V_nm * real(termo_acopl_p1) - real(termo_acopl_m1)
            soma_Fx += Q_nm * (Im_term_fx * D1 - Re_term_fx * D2)
            Re_term_fy = V_nm * real(termo_acopl_p1) + real(termo_acopl_m1)
            Im_term_fy = V_nm * imag(termo_acopl_p1) + imag(termo_acopl_m1)
            soma_Fy += Q_nm * (-Re_term_fy * D1 - Im_term_fy * D2)
        end
    end
    pre_fator = -(π * p.ρ₀ * p.ϕ₀_amp^2) / 2.0
    return pre_fator * soma_Fx, pre_fator * soma_Fy, pre_fator * soma_Fz
end