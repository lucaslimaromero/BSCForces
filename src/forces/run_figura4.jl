# run_figura4.jl

using Base.Threads
include("calculations.jl")
using Plots, LaTeXStrings


function gerar_grafico_figura4()
    println("Iniciando simulação para Figura 4 de Baresch et al. (2013)...")

    # --- 1. Parâmetros Fixos ---
    P₀ = 0.1e6    # 0.1 MPa
    c₀ = 1500.0   # m/s
    ρ₀ = 1000.0   # kg/m³
    cL_s = 2350.0 # m/s
    cT_s = 1120.0 # m/s
    ρ_s = 1080.0  # kg/m³
    f = 1e6       # 1 MHz
    ν = 1         # Carga topológica (helicity)
    
    n_max = 40 # Ordem máxima da expansão, suficiente para convergência
    a_div_λ_range = 0.01:0.0005:0.5
    #a_div_λ_range = 0.01:0.005:0.5
    α_axicon_lista_graus = [50.0, 60.0, 70.0]

    # --- 2. Preparar o Gráfico ---
    plot_final = plot(
        xlabel=L"a/\lambda",
        ylabel=L"F_z [\mathrm{N}]",
        title="Força Axial em Esfera de Poliestireno (ν=1, no eixo)",
        legend=:topleft,
        framestyle=:box,
        xlims=(0, 0.5),
        fontfamily="Computer Modern"
    )

    # --- 3. Loop Principal ---
    for α_graus in α_axicon_lista_graus
        println("Calculando para α = $(α_graus)°...")
        forcas_axiais = zeros(Float64, length(a_div_λ_range))
        
        # Pré-alocar matriz de BSCs
        bscs_2D = zeros(ComplexF64, n_max + 2, 2 * n_max + 3)

        for (i, a_div_λ) in enumerate(a_div_λ_range)
            # Cria a struct de parâmetros para esta iteração
            p = ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_graus, a_div_λ)

            # Calcula os coeficientes de espalhamento (dependem de 'a')
            α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, p)

            # Preenche a matriz de BSCs. Para o caso no eixo, é esparsa.
            # Apenas m=ν é não-nulo.
            fill!(bscs_2D, 0.0 + 0.0im) # Limpa a matriz
            for n in 0:n_max+1
                 # O índice de m é deslocado para ser positivo: m_idx = m + n_max + 2
                bscs_2D[n + 1, p.ν + n_max + 2] = bsccalc(n, p.ν, p, ρ_desloc=0.0)
            end
            
            # Calcula as três componentes da força
            ~, ~, Fz = calcular_forcas(p, n_max, bscs_2D, α_coeffs, β_coeffs)
            forcas_axiais[i] = Fz
            
            print("\rProgresso para α=$(α_graus)°: $(round(i/length(a_div_λ_range)*100, digits=1))%")
        end
        
        # Adiciona a curva ao gráfico
        plot!(plot_final, a_div_λ_range, forcas_axiais, 
            label="β = $(Int(α_graus))°", # β no artigo é o nosso α
            lw=2.5, linestyle=(α_graus==50 ? :solid : (α_graus==60 ? :dash : :dot))
        )
    end

    # --- 4. Finalizar e Salvar ---
    println("\nCálculo concluído. Gerando o gráfico...")
    savefig(plot_final, "figura4_baresch_reproduzida.png")
    println("Gráfico salvo como 'figura4_baresch_reproduzida.png'")
end

function encontrar_resonancia_beta_50()
    println("Iniciando busca de ressonância para β = 50°...")

    # Parâmetros fixos
    P₀ = 0.1e6; c₀ = 1500.0; ρ₀ = 1000.0;
    cL_s = 2350.0; cT_s = 1120.0; ρ_s = 1080.0;
    f = 1e6; ν = 1; n_max = 40
    a_div_λ_range = 0.01:0.0005:0.5 # Mais pontos para precisão

    α_graus = 50.0 # Foco apenas em β = 50°
    println("Calculando para β = $(α_graus)°...")

    # Inicializar variáveis para encontrar o máximo
    max_Fz = -Inf
    max_a_div_λ = 0.0

    forcas_axiais = zeros(Float64, length(a_div_λ_range))

    for (i, a_div_λ) in enumerate(a_div_λ_range)
        # Cria a struct de parâmetros para esta iteração
        p = ParametrosSimulacao(f, P₀, c₀, ρ₀, cL_s, cT_s, ρ_s, ν, α_graus, a_div_λ)

        # Calcula os coeficientes de espalhamento
        α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, p)

        # Preenche a matriz de BSCs
        bscs_2D = zeros(ComplexF64, n_max + 2, 2 * n_max + 3)
        fill!(bscs_2D, 0.0 + 0.0im)
        for n in 0:n_max+1
            bscs_2D[n + 1, p.ν + n_max + 2] = bsccalc(n, p.ν, p, ρ_desloc=0.0)
        end

        # Calcula a força axial
        ~, ~, Fz = calcular_forcas(p, n_max, bscs_2D, α_coeffs, β_coeffs)
        forcas_axiais[i] = Fz

        # Atualiza o máximo
        if Fz > max_Fz
            max_Fz = Fz
            max_a_div_λ = a_div_λ
        end
    end

    println("\nMáximo encontrado para β = 50°:")
    println("a/λ = $(max_a_div_λ), F_z = $(max_Fz) N")
end

# Executar a simulação para a Figura 4
gerar_grafico_figura4()

# Executar a busca de ressonância
encontrar_resonancia_beta_50()