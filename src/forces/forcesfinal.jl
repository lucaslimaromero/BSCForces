# --- Bibliotecas Necessárias ---
using SpecialFunctions
using LegendrePolynomials
using Plots
using LaTeXStrings

# ==============================================================================
# PARTE 1: ESTRUTURA PARA ORGANIZAR OS PARÂMETROS
# ==============================================================================

"""
Uma struct para agrupar todos os parâmetros físicos da simulação.
Isso torna o código mais limpo e fácil de gerenciar.
"""
struct ParametrosSimulacao
    # Meio (Água)
    f::Float64      # Frequência
    c₀::Float64     # Velocidade do som no meio
    ρ₀::Float64     # Densidade do meio
    P₀::Float64     # Amplitude de Pressão
    # Partícula (Poliestireno)
    a::Float64      # Raio da esfera
    ρ_s::Float64    # Densidade da esfera
    cL_s::Float64   # Velocidade longitudinal na esfera
    cT_s::Float64   # Velocidade transversal na esfera
    # Feixe
    ν::Int          # Carga topológica (ordem do feixe)
    α_axicon::Float64 # Ângulo do cone
    # Derivados (calculados automaticamente)
    λ::Float64      # Comprimento de onda
    k::Float64      # Número de onda
    ω::Float64      # Frequência angular
    ϕ₀_amp::Float64 # Amplitude do potencial de velocidade
end

# "Construtor" para facilitar a criação da struct
function ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_axicon_deg, a_div_λ)
    λ = c₀ / f
    a = a_div_λ * λ
    k = 2π / λ
    ω = 2π * f
    ϕ₀_amp = P₀ / (ω * ρ₀)
    α_axicon_rad = deg2rad(α_axicon_deg)
    
    return ParametrosSimulacao(f, c₀, ρ₀, P₀, a, ρ_s, cL_s, cT_s, ν, α_axicon_rad, λ, k, ω, ϕ₀_amp)
end


# ==============================================================================
# PARTE 2: SUAS FUNÇÕES DE CÁLCULO (MODIFICADAS PARA USAR A STRUCT)
# ==============================================================================

# Sua função, agora recebe a struct 'p'
function calcular_coeficientes_espalhamento(n_max::Int, p::ParametrosSimulacao)
    ω, k₀, a = p.ω, p.k, p.a
    kL = ω / p.cL_s
    kT = ω / p.cT_s
    
    x = k₀ * a
    xL = kL * a
    xT = kT * a

    α_coeffs = zeros(Float64, n_max + 1)
    β_coeffs = zeros(Float64, n_max + 1)

    sphericalhankelh(n, z) = sphericalbesselj(n, z) + im * sphericalbessely(n, z)
    ξ(n, z, h_func) = Float64(n) - z * h_func(n+1, z) / h_func(n, z)

    for n in 0:n_max
        ξ_n_x_j = ξ(n, x, sphericalbesselj)
        ξ_n_x_h = ξ(n, x, sphericalhankelh)

        dens_ratio = p.ρ_s / p.ρ₀
        
        ξ_n_xL_j = ξ(n, xL, sphericalbesselj)
        ξ_n_xT_j = ξ(n, xT, sphericalbesselj)

        termo1_num = 2 * n * (n + 1)
        termo1_den = (n^2 + n) - 0.5*xT^2 + 2*(n^2-1) / ( (ξ_n_xT_j - 1)/(ξ_n_xL_j - 1) - 1 )
        
        numerador_Rn = ξ_n_x_j * ( dens_ratio * ξ_n_xL_j + termo1_num / termo1_den ) - n*(n+1)*dens_ratio
        denominador_Rn = ξ_n_x_h * ( dens_ratio * ξ_n_xL_j + termo1_num / termo1_den ) - n*(n+1)*dens_ratio
        
        Rn = abs(denominador_Rn) < 1e-30 ? 0.0 + 0.0im : -numerador_Rn / denominador_Rn
        
        α_coeffs[n+1] = real(Rn)
        β_coeffs[n+1] = imag(Rn)
    end
    
    return α_coeffs, β_coeffs
end

# Sua função, agora recebe a struct 'p'
function bsccalc_deslocado(n::Int, m::Int, p::ParametrosSimulacao, ρ₀::Float64, ϕ₀::Float64)
    if n < abs(m); return ComplexF64(0.0); end
    
    kρ = p.k * sin(p.α_axicon)
    m_abs = abs(m)
    
    log_ratio_factorials = lgamma(n - m_abs + 1) - lgamma(n + m_abs + 1)
    on_axis_part = im^(n - m) * (2n + 1) * exp(log_ratio_factorials) * Plm(cos(p.α_axicon), n, m_abs)
    
    bessel_desloc = besselj(m - p.ν, kρ * ρ₀)
    fase_desloc_phi = cis(-(m - p.ν) * ϕ₀)
    
    return on_axis_part * bessel_desloc * fase_desloc_phi
end

# Sua função, agora recebe a struct 'p'
function calcular_forcas_transversais(p::ParametrosSimulacao, n_max, bscs_2D, α_coeffs, β_coeffs)
    soma_Fx = 0.0
    soma_Fy = 0.0

    for n in 0:n_max-1
        αn, βn = α_coeffs[n+1], β_coeffs[n+1]
        αn_p1, βn_p1 = α_coeffs[n+2], β_coeffs[n+2]

        D1 = αn + αn_p1 + 2 * (αn * αn_p1 + βn * βn_p1)
        D2 = βn_p1 - βn + 2 * (βn_p1 * αn - αn_p1 * βn)

        for m in -n:n
            Q_nm = 2 * exp(lgamma(n + m + 1) - lgamma(n - m + 1) - log(2n + 1) - log(2n + 3))
            V_nm = Float64((n+m+1)*(n+m+2))

            Anm = bscs_2D[n + 1, m + n_max + 1]
            Anp1_mp1 = (m + 1 <= n + 1) ? bscs_2D[n + 2, (m + 1) + n_max + 1] : 0.0
            Anp1_mm1 = (m - 1 >= -n - 1) ? bscs_2D[n + 2, (m - 1) + n_max + 1] : 0.0

            Im_term_fx = V_nm * imag(Anm * conj(Anp1_mp1)) - imag(Anm * conj(Anp1_mm1))
            Re_term_fx = V_nm * real(Anm * conj(Anp1_mp1)) - real(Anm * conj(Anp1_mm1))
            soma_Fx += Q_nm * (Im_term_fx * D1 - Re_term_fx * D2)
            
            Re_term_fy = V_nm * real(Anm * conj(Anp1_mp1)) + real(Anm * conj(Anp1_mm1))
            Im_term_fy = V_nm * imag(Anm * conj(Anp1_mp1)) + imag(Anm * conj(Anp1_mm1))
            soma_Fy += Q_nm * (-Re_term_fy * D1 - Im_term_fy * D2)
        end
    end

    pre_fator = (π * p.ρ₀ * p.ϕ₀_amp^2) / 2.0
    return pre_fator * soma_Fx, pre_fator * soma_Fy
end


# ==============================================================================
# PARTE 3: SCRIPT PRINCIPAL PARA GERAR O GRÁFICO 
# ==============================================================================

function gerar_grafico_figura2()
    println("Iniciando simulação para Figura 2 de Baresch et al. (2013)...")

    # --- 1. Definindo Parâmetros ---
    p = ParametrosSimulacao(
        1e6,      # f (Hz)
        0.1e6,    # P₀ (Pa)
        1500.0,   # c₀ (m/s)
        1000.0,   # ρ₀ (kg/m³)
        2350.0,   # cL_s (m/s)
        1120.0,   # cT_s (m/s)
        1080.0,   # ρ_s (kg/m³)
        1,        # ν (ordem do feixe)
        25.0,     # α_axicon (graus)
        0.01      # a/λ
    )
    
    n_max = 50
    desloc_div_λ = 0.0:0.05:4.0
    deslocamentos_ρ = desloc_div_λ * p.λ

    # --- 2. Pré-cálculo dos Coeficientes de Espalhamento (feito uma vez) ---
    println("Calculando coeficientes de espalhamento...")
    α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, p)

    # --- 3. Loop de Simulação Otimizado ---
    forcas_radiais = zeros(Float64, length(deslocamentos_ρ))
    
    # OTIMIZAÇÃO: Criamos a matriz de BSCs UMA VEZ, fora do loop.
    bscs_2D = zeros(ComplexF64, n_max + 2, 2 * n_max + 1)
    
    println("Iniciando loop de cálculo da força...")
    for (i, ρ₀_desloc) in enumerate(deslocamentos_ρ)
        ϕ₀_desloc = 0.0

        # ATUALIZAÇÃO: Em vez de criar, apenas preenchemos a matriz já existente.
        for n in 0:n_max+1
            for m in -n:n
                if abs(m) > n_max; continue; end
                bscs_2D[n + 1, m + n_max + 1] = bsccalc_deslocado(n, m, p, ρ₀_desloc, ϕ₀_desloc)
            end
        end

        Fx, Fy = calcular_forcas_transversais(p, n_max, bscs_2D, α_coeffs, β_coeffs)
        forcas_radiais[i] = Fx
        
        print("\rProgresso: $(round(i/length(deslocamentos_ρ)*100, digits=1))%")
    end

    # --- 4. Plotagem (sem alterações) ---
    println("\nCálculo concluído. Gerando o gráfico...")
    
    plot(desloc_div_λ, forcas_radiais .* 1e12,
        label="Cálculo (Ondas Parciais)",
        lw=2,
        xlabel="ρ/λ",
        ylabel="F_ρ [pN]",
        title="Força Radial em Esfera de Poliestireno (ν=1, a=0.01λ)",
        legend=:bottomright,
        framestyle=:box
    )
    
    perfil_bessel = (x -> abs2(besselj(p.ν, p.k * sin(p.α_axicon) * x * p.λ)))
    plot!(desloc_div_λ, x -> 2.0 * perfil_bessel(x),
        label="Perfil |J₁|² (arbitrário)",
        linestyle=:dot,
        color=:gray,
        lw=1.5
    )

    savefig("figura2_baresch_refatorada.png")
    println("Gráfico salvo como 'figura2_baresch_refatorada.png'")
end

# Executar a simulação
gerar_grafico_figura2()