# run_figura5.jl

include("calculations.jl")
using Plots, LaTeXStrings

function gerar_grafico_figura5()
    println("Iniciando simulação para Figura 5 de Baresch et al. (2013)...")

    # --- 1. Definir os dois conjuntos de parâmetros para as curvas da Figura 5 ---
    # Parâmetros extraídos diretamente da legenda da Fig. 5 do artigo.
    params_caso1 = (β = 50.0, a_div_λ = 0.30, estilo = :solid, cor = :black, label = "a/λ=0.300 (Não ressonante)")
    params_caso2 = (β = 50.0, a_div_λ = 0.358, estilo = :dash, cor = :black, label = "a/λ=0.358 (Ressonante)")
    
    casos = [params_caso1, params_caso2]
    
    # --- 2. Preparar o gráfico ---
    plot_final = plot(
        xlabel=L"\rho/\lambda", 
        ylabel=L"F_\rho [N]", 
        title="Força Radial Transversal (Figura 5)",
        legend=:bottomright,
        fontfamily="Computer Modern",
        framestyle=:box
    )

    # --- 3. Loop principal para calcular as curvas ---
    for caso in casos
        println("Calculando para: $(caso.label)...")
        
        p = ParametrosSimulacao(1e6, 0.1e6, 1500.0, 1000.0, 2350.0, 1120.0, 1080.0, 1, caso.β, caso.a_div_λ)
        
        # Para estes tamanhos de partícula, DEVEMOS usar a função de cálculo geral.
        n_max = 40 # Usamos um n_max maior para garantir a convergência
        
        desloc_div_λ = 0.0:0.02:1.5
        #desloc_div_λ = 0.0:0.004:1.5
        deslocamentos_ρ = desloc_div_λ * p.λ
        forcas_radiais = zeros(Float64, length(deslocamentos_ρ))

        # Pré-cálculo dos coeficientes de espalhamento (feito uma vez por caso)
        # Chamando a função correta para este regime!
        α_coeffs, β_coeffs = calcular_coeficientes_espalhamento(n_max, p)
        
        bscs_2D = zeros(ComplexF64, n_max + 2, 2 * n_max + 3)

        for (i, ρ₀_desloc) in enumerate(deslocamentos_ρ)
            # Preenche a matriz de BSCs para este deslocamento
            for n in 0:n_max + 1, m in -n:n
                bscs_2D[n + 1, m + n_max + 2] = bsccalc(n, m, p, ρ_desloc=ρ₀_desloc)
            end

            # Calcula as forças usando a função unificada, mas só queremos Fx nesse plot
            Fx, ~, ~ = calcular_forcas(p, n_max, bscs_2D, α_coeffs, β_coeffs)
            
            # Para ajustar o referencial para o feixe parado com a partícula "em movimento", invertemos a força
            forcas_radiais[i] = -Fx

            print("\rProgresso para $(caso.label): $(round(i/length(deslocamentos_ρ)*100, digits=1))%")
        end
        println()

        # Adiciona a curva ao gráfico. Multiplicamos por 1e7 para a escala do artigo.
        plot!(plot_final, desloc_div_λ, forcas_radiais .* 1e7,
            label=caso.label, lw=2.5, linestyle=caso.estilo, color=caso.cor)
    end

    # --- 4. Finalizar e salvar o gráfico ---
    # Adicionando a linha do zero e ajustando os eixos para bater com a referência
    hline!(plot_final, [0], color=:black, lw=1, label="")
    #ylims!(plot_final, (-2, 1)) # Limites originais
    #yticks!(plot_final, -2:0.5:1)
    
    # Adiciona o 'x 10⁻⁷' no eixo Y, como no artigo
    ylabel!(plot_final, L"F_\rho[\mathrm{N}] \quad (\,\times\,10^{-7})")
    
    savefig(plot_final, "figura5_baresch_reproduzida.png")
    println("\nGráfico final para a Figura 5 salvo como 'figura5_baresch_reproduzida.png'")
end

gerar_grafico_figura5()