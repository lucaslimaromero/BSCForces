using SpecialFunctions, LegendrePolynomials, Plots, Base.Threads, LaTeXStrings, Plots.PlotMeasures

# Sua função bsccalc original (sem alterações aqui)
function bsccalc(n::Int, m::Int, α::Float64) 
    if n < abs(m) 
        return ComplexF64(0.0)
    end
    m_abs = abs(m)
    log_fact_nm_plus_1 = lgamma(Float64(n - m_abs + 1))
    log_fact_npm_plus_1 = lgamma(Float64(n + m_abs + 1))
    ratio_factorials = exp(log_fact_nm_plus_1 - log_fact_npm_plus_1)
    phase_factor = ComplexF64(im^(n - m))
    amplitude_factor = (2.0 * n + 1.0)
    legendre_val = Plm(cos(α), n, m_abs) # Plm da biblioteca LegendrePolynomials
    return phase_factor * amplitude_factor * ratio_factorials * legendre_val
end

# Função otimizada para calcular a amplitude da onda parcial
function optimized_partialwave_calculator(
    psiamp_val::Number, 
    order_val::Int, 
    r_coord::Float64, 
    k_const::Float64, 
    phi_coord::Float64, 
    n_max_val::Int, 
    precomputed_bscs_arr::Vector{ComplexF64}, # BSCs para o order_val e axiconang atuais
    precomputed_plm_at_pi_2_arr::Vector{Float64} # Plm(0, n, abs(order_val))
)
    # Tratamento do ponto central r=0
    if r_coord == 0.0
        return order_val == 0 ? ComplexF64(psiamp_val) : ComplexF64(0.0)
    end

    kr_val = k_const * r_coord
    psi_total = ComplexF64(0.0)

    # Loop principal apenas sobre n, usando valores pré-calculados
    # O loop deve ir de n = abs(order_val) até n_max_val
    # Pois Plm(0, n, abs(order_val)) requer n >= abs(order_val)
    # E bsccalc também tem a condição n >= abs(m)
    for n_idx in abs(order_val):n_max_val
        # Acessa tabelas pré-calculadas (ajuste de índice +1 pois n_idx pode ser 0)
        term_bsc = precomputed_bscs_arr[n_idx + 1]
        term_plm = precomputed_plm_at_pi_2_arr[n_idx + 1]

        # Se Plm(0,n,m) for zero (porque n+m é ímpar), term_plm será zero, anulando a contribuição
        if term_plm == 0.0 # Otimização adicional: pular cálculo se Plm é zero
            continue
        end
        if term_bsc == ComplexF64(0.0) # Se BSC for zero (ex: n < abs(order_val))
             continue
        end


        spher_val = SpecialFunctions.sphericalbesselj(n_idx, kr_val) # kr_val depende de r_coord (pixel)
        
        psi_total += term_bsc * spher_val * term_plm * cis(order_val * phi_coord)
    end
    return psi_total * psiamp_val
end


function makeplotpartial_optimized(nx::Int, ny::Int, rx_init::Float64, ry_init::Float64, ψ₀_amp::Number)
    # z = 0 é implícito e usado para θ = π/2
    println("Threads disponíveis: ", Threads.nthreads())
    k_val = 1000.0 * 2.0 * pi / 1.54
    n_max = 1000

    axicon_angles_rad = [deg2rad(1.0), deg2rad(10.0), deg2rad(40.0)]
    orders_list = [0, 1, 5]
    
    num_plots = length(axicon_angles_rad) * length(orders_list)
    heatmaps_collection = Vector{Any}(undef, num_plots)
    plot_labels = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ] # Garanta que tem rótulos suficientes
    current_plot_idx = 0

    # Variáveis para o reescalonamento de coordenadas (conforme seu código original)
    current_rx = rx_init
    current_ry = ry_init
    mag_factor = 10.0 # mag inicial do seu código original

    for alpha_rad_current in axicon_angles_rad
        # Recalcular coordenadas x e y se elas mudam por alpha_rad_current devido ao 'mag_factor'
        # No seu código original, x e y são reescalados *após* o loop de alpha_rad_current
        # Vou manter essa lógica, mas é importante estar ciente
        x_coords_vec = range(-current_rx, current_rx, length=nx)
        y_coords_vec = range(-current_ry, current_ry, length=ny)

        for order_current in orders_list
            current_plot_idx += 1
            bessel_beam_intensity_map = Array{Float64}(undef, nx, ny)

            # --- Início dos Pré-cálculos para (alpha_rad_current, order_current) ---
            # 1. Pré-calcular BSCs: bsccalc(n, order_current, alpha_rad_current)
            #    bsccalc internamente usa Plm(cos(alpha_rad_current), n, abs(order_current))
            precomputed_bscs_arr = Vector{ComplexF64}(undef, n_max + 1) # Para n = 0 to n_max
            for n_val in 0:n_max
                precomputed_bscs_arr[n_val + 1] = bsccalc(n_val, order_current, alpha_rad_current)
            end

            # 2. Pré-calcular Plm(0.0, n, abs(order_current)) pois cos(θ)=0 para z=0
            precomputed_plm_at_pi_2_arr = Vector{Float64}(undef, n_max + 1) # Para n = 0 to n_max
            abs_order_current = abs(order_current)
            for n_val in 0:n_max
                if n_val < abs_order_current || (n_val + abs_order_current) % 2 != 0 
                    # P_n^m(0) é zero se n < m ou se n+m é ímpar
                    precomputed_plm_at_pi_2_arr[n_val + 1] = 0.0
                else
                    precomputed_plm_at_pi_2_arr[n_val + 1] = Plm(0.0, n_val, abs_order_current)
                end
            end
            # --- Fim dos Pré-cálculos ---
            
            @threads for j_idx in 1:ny # Paralelização por linhas da imagem
                yi_coord = y_coords_vec[j_idx]
                for i_idx in 1:nx
                    xi_coord = x_coords_vec[i_idx]

                    rho_coord = sqrt(xi_coord^2 + yi_coord^2)
                    # r_coord é igual a rho_coord porque z=0

                    phi_coord = atan(yi_coord, xi_coord) # atan(y,x) é mais robusto
                    # Se precisar de phi em [0, 2π): if phi_coord < 0 phi_coord += 2*pi; end

                    amplitude_complexa = optimized_partialwave_calculator(
                        ψ₀_amp, order_current, rho_coord, k_val, phi_coord, 
                        n_max, precomputed_bscs_arr, precomputed_plm_at_pi_2_arr
                    )
                    # Usar abs2 para intensidade (equivale a abs(val)^2 mas pode ser mais rápido)
                    bessel_beam_intensity_map[i_idx, j_idx] = abs2(amplitude_complexa) 
                end
            end
            
            # Adiciona o heatmap ao vetor (usando x_coords_vec, y_coords_vec para os eixos)
            # Multiplicar por 1000 para mm se as coordenadas rx, ry estão em metros
            heatmap_title = latexstring("(", plot_labels[current_plot_idx], ") \\alpha=", round(rad2deg(alpha_rad_current),digits=1),"^\\circ, \\nu=", order_current)
            heatmaps_collection[current_plot_idx] = heatmap(1000 * x_coords_vec, 1000 * y_coords_vec, bessel_beam_intensity_map, 
                                                            xlabel=latexstring("x (mm)"), ylabel=latexstring("y (mm)"), 
                                                            title=heatmap_title, aspect_ratio=:equal,
                                                            left_margin=5mm, bottom_margin=5mm) # Ajuste de margens
        end # fim loop order_current
        
        # Lógica de reescalonamento do seu código original para a *próxima* iteração de alpha_rad_current
        current_rx = current_rx / mag_factor
        current_ry = current_ry / mag_factor
        mag_factor = 4.0 # mag para as próximas iterações
    end # fim loop alpha_rad_current
    
    # Plotar todos os heatmaps juntos
    final_plot = plot(heatmaps_collection..., layout=(length(axicon_angles_rad), length(orders_list)), size=(1200, 900), titlefontsize=10)
    savefig(final_plot, "partial_wave_exp1_optimized.png")
    println("Simulação otimizada concluída. Imagem salva em partial_wave_exp1_optimized.png")
end

# Para executar:
# Certifique-se de que o Julia foi iniciado com threads, ex: julia -t auto ou julia -t 4
@time makeplotpartial_optimized(100, 100, 0.2, 0.2, 1.0)
# Para testes mais rápidos, reduza nx, ny ou n_max:
# @time makeplotpartial_optimized(10, 10, 0.2, 0.2, 1.0) # Resolução bem baixa para teste rápido
# Para a sua resolução original e teste de tempo:
# @time makeplotpartial_optimized(30, 30, 0.2, 0.2, 1.0) # nx=30, ny=30
# Se quiser maior resolução para a imagem final (vai demorar mais):
# @time makeplotpartial_optimized(100, 100, 0.2, 0.2, 1.0) # nx=100, ny=100