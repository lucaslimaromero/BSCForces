# --- Bibliotecas Necessárias ---
using SpecialFunctions
using LegendrePolynomials
using Plots
using LaTeXStrings
using Plots.PlotMeasures

# ==============================================================================
# PARTE 1: CÁLCULO DOS COEFICIENTES DE ESPALHAMENTO (VERSÃO FINAL E CORRIGIDA)
# Esta versão é baseada na formulação clássica de Faran/Anderson e é mais estável.
# ==============================================================================
function calcular_coeficientes_espalhamento(n_max::Int, a::Float64, f::Float64, ρ₀::Float64, c₀::Float64, ρ_s::Float64, cL_s::Float64, cT_s::Float64)
    ω = 2π * f
    k₀ = ω / c₀
    kL = ω / cL_s
    kT = ω / cT_s
    
    x = k₀ * a
    xL = kL * a
    xT = kT * a

    α_coeffs = zeros(Float64, n_max + 1)
    β_coeffs = zeros(Float64, n_max + 1)

    # Funções de Hankel
    sphericalhankelh(n, z) = sphericalbesselj(n, z) + im * sphericalbessely(n, z)

    # Função auxiliar para a razão ξ_n(z) = z * h_n'(z) / h_n(z)
    # Usa a relação de recorrência para evitar calcular a derivada diretamente
    ξ(n, z, h_func) = Float64(n) - z * h_func(n+1, z) / h_func(n, z)

    for n in 0:n_max
        # Termos que dependem do fluido (água)
        ξ_n_x_j = ξ(n, x, sphericalbesselj)
        ξ_n_x_h = ξ(n, x, sphericalhankelh)

        # Termos que dependem do sólido (poliestireno)
        dens_ratio = ρ_s / ρ₀
        
        ξ_n_xL_j = ξ(n, xL, sphericalbesselj)
        ξ_n_xT_j = ξ(n, xT, sphericalbesselj)

        # Fórmula de Anderson (1950)
        termo1_num = 2 * n * (n + 1)
        termo1_den = (n^2 + n) - 0.5*xT^2 + 2*(n^2-1) / ( (ξ_n_xT_j - 1)/(ξ_n_xL_j - 1) - 1 )
        
        numerador_Rn = ξ_n_x_j * ( dens_ratio * ξ_n_xL_j + termo1_num / termo1_den ) - n*(n+1)*dens_ratio
        denominador_Rn = ξ_n_x_h * ( dens_ratio * ξ_n_xL_j + termo1_num / termo1_den ) - n*(n+1)*dens_ratio
        
        # O coeficiente de espalhamento é R_n = -numerador / denominador
        if abs(denominador_Rn) < 1e-30
            Rn = 0.0 + 0.0im
        else
            Rn = -numerador_Rn / denominador_Rn
        end
        
        α_coeffs[n+1] = real(Rn)
        β_coeffs[n+1] = imag(Rn)
    end
    
    return α_coeffs, β_coeffs
end

# ==============================================================================
# PARTE 1: CÁLCULO DOS COEFICIENTES DE ESPALHAMENTO (APROXIMAÇÃO DE RAYLEIGH)
# Para a << λ, esta é a forma mais simples e correta de calcular os coeficientes.
# ==============================================================================
function calcular_coeffs_rayleigh(n_max::Int, a::Float64, k::Float64, ρ₀::Float64, c₀::Float64, ρ_s::Float64, cL_s::Float64)
    α_coeffs = zeros(Float64, n_max + 1)
    β_coeffs = zeros(Float64, n_max + 1)
    
    x = k * a # ka
    
    # Coeficiente de compressibilidade (K = ρ * c²)
    K_fluido = ρ₀ * c₀^2
    K_solido = ρ_s * cL_s^2
    
    # Termo f₀ para o espalhamento monopolar (n=0)
    f₀ = 1.0 - (K_fluido / K_solido)
    
    # Termo f₁ para o espalhamento dipolar (n=1)
    f₁ = (2.0 * (ρ_s - ρ₀)) / (2.0 * ρ_s + ρ₀)
    
    # --- Coeficiente n=0 ---
    # R₀ ≈ -i/3 * (ka)³ * f₀
    if n_max >= 0
        R0 = -im * (x^3 / 3.0) * f₀
        α_coeffs[1] = real(R0)
        β_coeffs[1] = imag(R0)
    end

    # --- Coeficiente n=1 ---
    # R₁ ≈ i/3 * (ka)³ * f₁  (Nota: há um sinal de - que cancela com outro na derivação completa)
    # A convenção mais comum para a força resulta em um sinal positivo aqui.
    if n_max >= 1
        #R1 = -im * (x^3 / 3.0) * f₁
        R1 = +im * (x^3 / 3.0) * f₁
        α_coeffs[2] = real(R1)
        β_coeffs[2] = imag(R1)
    end

    # Para n > 1, os coeficientes são considerados zero nesta aproximação.
    
    return α_coeffs, β_coeffs
end

# ==============================================================================
# PARTE 2: CÁLCULO DOS BSCs PARA FEIXE DESLOCADO
# ==============================================================================

"""
Calcula UM Coeficiente de Forma do Feixe (BSC) para um feixe de Bessel de ordem ν,
que está DESLOCADO da origem.
"""
function bsccalc_deslocado(n::Int, m::Int, α_axicon::Float64, ν::Int, k::Float64, ρ₀::Float64, ϕ₀::Float64)
    if n < abs(m); return ComplexF64(0.0); end
    
    kρ = k * sin(α_axicon)
    m_abs = abs(m)
    
    log_ratio_factorials = lgamma(n - m_abs + 1) - lgamma(n + m_abs + 1)
    on_axis_part = ComplexF64(im)^(n - m) * (2n + 1) * exp(log_ratio_factorials) * Plm(cos(α_axicon), n, m_abs)
    
    bessel_desloc = besselj(m - ν, kρ * ρ₀)
    fase_desloc_phi = cis(-(m - ν) * ϕ₀)
    
    return on_axis_part * bessel_desloc * fase_desloc_phi
end


# ==============================================================================
# PARTE 3: IMPLEMENTAÇÃO DAS EQUAÇÕES DE FORÇA
# ==============================================================================

"""
Calcula as componentes cartesianas da força Fx e Fy, seguindo rigorosamente
a fórmula do projeto de IC (Eq. 6a e 6b).
"""
function calcular_forcas_transversais(ρ₀_fluid, ϕ₀_amp, n_max, bscs_2D, α_coeffs, β_coeffs)
    soma_Fx = ComplexF64(0.0)
    soma_Fy = ComplexF64(0.0)

    for n in 0:n_max-1
        αn, βn = α_coeffs[n+1], β_coeffs[n+1]
        αn_p1, βn_p1 = α_coeffs[n+2], β_coeffs[n+2]

        D1 = αn + αn_p1 + 2 * (αn * αn_p1 + βn * βn_p1)
        D2 = βn_p1 - βn + 2 * (βn_p1 * αn - αn_p1 * βn)

        # O limite da soma em m é |m|<=n, mas a fórmula usa termos A_{n+1}^{m+1},
        # então o loop em m precisa ser cuidadoso nos limites.
        for m in -(n+1):(n+1)
            # Pula iterações onde Qnm não seria válido.
            if abs(m) > n; continue; end

            m_idx = m + n_max + 1
            
            # Coeficientes Qnm e Vnm
            Q_nm = 2 * exp(lgamma(n + m + 1) - lgamma(n - m + 1) - log(2n + 1) - log(2n + 3))
            V_nm = Float64((n+m+1)*(n+m+2)) # Definição da Eq. 6 do seu projeto.

            # Acessa os BSCs necessários
            Anm = bscs_2D[n + 1, m_idx]
            
            # Para A_{n+1}^{m+1}
            Anp1_mp1 = (m + 1 <= n + 1) ? bscs_2D[n + 2, (m + 1) + n_max + 1] : 0.0
            
            # Para A_{n+1}^{m-1}
            Anp1_mm1 = (m - 1 >= -n - 1) ? bscs_2D[n + 2, (m - 1) + n_max + 1] : 0.0

            # --- Cálculo de Fx (Eq. 6a) ---
            Im_term1_fx = V_nm * imag(Anm * conj(Anp1_mp1)) - imag(Anm * conj(Anp1_mm1))
            Re_term1_fx = V_nm * real(Anm * conj(Anp1_mp1)) - real(Anm * conj(Anp1_mm1))
            soma_Fx += Q_nm * (Im_term1_fx * D1 - Re_term1_fx * D2)
            
            # --- Cálculo de Fy (Eq. 6b) ---
            Re_term2_fy = V_nm * real(Anm * conj(Anp1_mp1)) + real(Anm * conj(Anp1_mm1))
            Im_term2_fy = V_nm * imag(Anm * conj(Anp1_mp1)) + imag(Anm * conj(Anp1_mm1))
            soma_Fy += Q_nm * (-Re_term2_fy * D1 - Im_term2_fy * D2)
        end
    end

    #pre_fator = -(π * ρ₀_fluid * ϕ₀_amp^2) / 2.0
    pre_fator = (π * ρ₀_fluid * ϕ₀_amp^2) / 2.0
    return pre_fator * real(soma_Fx), pre_fator * real(soma_Fy)
end


# ==============================================================================
# PARTE 4: SCRIPT PRINCIPAL PARA GERAR O GRÁFICO
# ==============================================================================

function gerar_grafico_figura2()
    
    println("Iniciando simulação para Figura 2 de Baresch et al. (2013)...")

    # --- 1. Definir Parâmetros Físicos ---
    # [cite_start]Parâmetros da onda e do meio (água) [cite: 1]
    f = 1e6       # Frequência: 1 MHz
    P₀ = 0.1e6    # Amplitude de Pressão: 0.1 MPa
    c₀ = 1500.0   # Velocidade do som na água: 1500 m/s
    ρ₀ = 1000.0   # Densidade da água: 1000 kg/m^3
    
    # [cite_start]Parâmetros da esfera (poliestireno) [cite: 1]
    cL_s = 2350.0 # Velocidade longitudinal
    cT_s = 1120.0 # Velocidade transversal
    ρ_s = 1080.0  # Densidade
    
    # [cite_start]Parâmetros do feixe [cite: 1]
    ν = 1                     # Carga topológica (helicity m=1)
    α_axicon = deg2rad(25.0)  # Ângulo do cone (beta = 25 graus)

    # Parâmetros derivados
    λ = c₀ / f
    k = 2π / λ
    ω = 2π * f
    ϕ₀_amp = P₀ / (ω * ρ₀)
    a = 0.01 * λ  # Raio da esfera a = 0.01λ
    
    # Parâmetros de simulação
    n_max = 50  # Para a aproximação de Rayleigh, n pequenos dominam
    desloc_div_λ = 0.0:0.05:4.0 # Eixo x do gráfico (ρ/λ)
    deslocamentos_ρ = desloc_div_λ * λ

    # --- 2. Pré-cálculo dos Coeficientes de Espalhamento ---
    println("Calculando coeficientes de espalhamento...")
    
    # Linha original:
    α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, a, f, ρ₀, c₀, ρ_s, cL_s, cT_s)

    # --- 3. Loop de Simulação ---
    forcas_radiais = zeros(Float64, length(deslocamentos_ρ))
    println("Iniciando loop de cálculo da força...")

    for (i, ρ₀_desloc) in enumerate(deslocamentos_ρ)
        # Para simplificar, o deslocamento é sempre no eixo x (ϕ₀ = 0)
        ϕ₀_desloc = 0.0

        # Pré-calcula TODOS os BSCs para este deslocamento
        bscs_2D = zeros(ComplexF64, n_max + 1, 2 * n_max + 1)
        for n in 0:n_max
            for m in -n:n
                m_idx = m + n_max + 1
                bscs_2D[n + 1, m_idx] = bsccalc_deslocado(n, m, α_axicon, ν, k, ρ₀_desloc, ϕ₀_desloc)
            end
        end

        # Calcula as forças cartesianas Fx e Fy
        Fx, Fy = calcular_forcas_transversais(ρ₀, ϕ₀_amp, n_max, bscs_2D, α_coeffs, β_coeffs)
        
        # A força radial é F_ρ = Fx*cos(ϕ₀) + Fy*sin(ϕ₀). Como ϕ₀=0, F_ρ = Fx.
        forcas_radiais[i] = Fx
        
        print("\rProgresso: $(round(i/length(deslocamentos_ρ)*100, digits=1))%")
    end

    # --- 4. Plotagem ---
    println("\nCálculo concluído. Gerando o gráfico...")
    
    plot(desloc_div_λ, forcas_radiais .* 1e12,
        label="Cálculo (Ondas Parciais)",
        lw=2,
        xlabel="ρ/λ",
        ylabel="F_ρ [N]  (x 10⁻¹²)",
        title="Força Radial em Esfera de Poliestireno (ν=1, a=0.01λ)",
        legend=:bottomright,
        framestyle=:box
    )
    
    # Adicionar a curva tracejada do perfil do feixe (aproximação)
    perfil_bessel = (x -> abs2(besselj(ν, k * sin(α_axicon) * x * λ)))
    plot!(desloc_div_λ, x -> 2.0 * perfil_bessel(x),
        label="Perfil |J₁|² (arbitrário)",
        linestyle=:dot,
        color=:gray,
        lw=1.5
    )

    savefig("figura2_baresch_reproduzida.png")
    println("Gráfico salvo como 'figura2_baresch_reproduzida.png'")
end

# ==============================================================================
# FUNÇÃO DE TESTE: Coeficientes para uma Esfera Rígida e Fixa
# Esta fórmula é muito mais simples e serve para verificarmos o resto do código.
# ==============================================================================
function calcular_coeffs_esfera_rigida_TESTE(n_max::Int, a::Float64, k::Float64)
    
    
    α_coeffs = zeros(Float64, n_max + 1)
    β_coeffs = zeros(Float64, n_max + 1)
    x = k * a

    # Derivadas das funções de bessel esféricas
    sphericalbesselj_deriv(n, z) = sphericalbesselj(n-1, z) - (n+1)/z * sphericalbesselj(n, z)
    sphericalhankelh_deriv(n, z) = (sphericalbesselj(n-1, z) - (n+1)/z * sphericalbesselj(n, z)) + im * (sphericalbessely(n-1, z) - (n+1)/z * sphericalbessely(n, z))

    for n in 0:n_max
        # Para uma esfera rígida, Rn = -j_n'(x) / h_n'(x)
        numerador = sphericalbesselj_deriv(n, x)
        denominador = sphericalhankelh_deriv(n, x)

        if abs(denominador) < 1e-30
            Rn = 0.0 + 0.0im
        else
            Rn = -numerador / denominador
        end
        
        α_coeffs[n+1] = real(Rn)
        β_coeffs[n+1] = imag(Rn)
    end
    
    return α_coeffs, β_coeffs
end

# Executar a simulação
gerar_grafico_figura2()