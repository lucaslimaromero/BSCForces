# run_figura2.jl

include("calculations.jl")
using Plots, LaTeXStrings

function gerar_grafico_figura2()
    println("Iniciando simulação para Figura 2 de Baresch et al. (2013)...")

    # --- 1. Parâmetros ---
    p = ParametrosSimulacao(
        1e6,      # f (Hz)
        0.1e6,    # P₀ (Pa)
        1500.0,   # c₀ (m/s)
        1000.0,   # ρ₀ (kg/m³)
        2350.0,   # cL_s (m/s)
        1120.0,   # cT_s (m/s)
        1080.0,   # ρ_s (kg/m³)
        1,        # ν (ordem do feixe)
        25.0,     # α_axicon (graus, β no artigo)
        0.01      # a/λ (fixo)
    )
    
    # =================================================================
    # A CORREÇÃO ESTÁ AQUI:
    # Para partículas pequenas (regime de Rayleigh), apenas os primeiros
    # termos da série são significativos. Usar um n_max grande causa
    # instabilidade numérica (divisão por zero).
    n_max = 5 
    # =================================================================
    
    desloc_div_λ = 0.0:0.05:4.0
    deslocamentos_ρ = desloc_div_λ * p.λ

    # --- 2. Cálculos ---
    println("Calculando coeficientes de espalhamento (uma vez)...")
    #α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, p)
    α_coeffs, β_coeffs = calcular_coeffs_rayleigh(n_max, p)
    
    forcas_radiais = zeros(Float64, length(deslocamentos_ρ))
    bscs_2D = zeros(ComplexF64, n_max + 2, 2 * n_max + 3)
    
    println("Iniciando loop de cálculo da força...")
    for (i, ρ₀_desloc) in enumerate(deslocamentos_ρ)
        ϕ₀_desloc = 0.0 # Deslocamento apenas no eixo x

        for n in 0:n_max+1
            for m in -n:n
                bscs_2D[n + 1, m + n_max + 2] = bsccalc(n, m, p, ρ_desloc=ρ₀_desloc, ϕ_desloc=ϕ₀_desloc)
            end
        end

        Fx, ~, ~ = calcular_forcas(p, n_max, bscs_2D, α_coeffs, β_coeffs)
        forcas_radiais[i] = Fx 
        
        print("\rProgresso: $(round(i/length(deslocamentos_ρ)*100, digits=1))%")
    end

    # --- 3. Plotagem ---
    println("\nCálculo concluído. Gerando o gráfico...")
    
    plot(desloc_div_λ, forcas_radiais .* 1e12, # Convertendo para piconewtons (pN)
        label="Força Calculada",
        lw=2.5,
        color=:black,
        seriestype=:line, # Apenas a linha
        xlabel=L"\rho/\lambda",
        ylabel=L"F_\rho[pN]",
        title="Força Radial (a/λ=0.01, β=25°)",
        legend=:bottomright,
        framestyle=:box,
        fontfamily="Computer Modern"
    )

    # Perfil de intensidade para comparação
    k_rho = p.k * sin(p.α_axicon)
    perfil_bessel_int = (rho -> abs2(besselj(p.ν, k_rho * rho)))
    max_forca_abs = maximum(abs, forcas_radiais)
    # A escala arbitrária é só para o perfil de intensidade caber bem no gráfico
    plot!(desloc_div_λ, x -> 1.5e-12 * perfil_bessel_int.(x * p.λ) * 1e12,
        label="Perfil "*L"|J_1|^2",
        linestyle=:dot,
        color=:gray,
        lw=2
    )

    hline!([0], color=:black, lw=1, label="")
    savefig("figura2_baresch_reproduzida.png")
    println("Gráfico salvo como 'figura2_baresch_reproduzida.png'")
end

gerar_grafico_figura2()