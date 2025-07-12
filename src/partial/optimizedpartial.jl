#using Bessels # Alta precisão para funções de Bessel
using SpecialFunctions # inclui funções como besselj e sphericalbesselj
using LegendrePolynomials
using Plots
using BenchmarkTools
using Test 
using LaTeXStrings
using Plots.PlotMeasures
using Base.Threads

# α: axiconang
# h = kρ = k * sinα      
# β = kz = k * cosα      

# A_{n,m} = i^{(n - m)} (2n+1) * [(n-|m|)! / (n+|m|)!] * P_n^{|m|}(cos(α))
function bsccalc(n::Int, m::Int, α::Float64) 
    if n < abs(m) 
        return ComplexF64(0.0)
    end
    # Módulo de m
    m_abs = abs(m)

    # Γ(n+1) = n! 
    # Aplicando log dos dois lados: ln(Γ(n+1)) = ln(n!)
    #  = lgamma(n+1)

    # Precisamos de (n-m)!/(n+m)!

    # Isso, pela nossa função gamma seria (n-m)!/(n+m)! = Γ(n-m+1)/Γ(n+m+1)
    # Aplicando ln dos dois lados
    # ln[(n-m)!/(n+m)!] = ln(Γ(n-m+1)/Γ(n+m+1)). Pela propriedade de logarítmos ln(a/b) = ln(a) - ln(b)

    num = lgamma(Float64(n - m_abs + 1)) # ln(Γ(n-m+1)
    denom = lgamma(Float64(n + m_abs + 1)) # ln(Γ(n+m+1))

    # ln[(n-m)!/(n+m)!] = ln(Γ(n-m+1) - ln(Γ(n+m+1))
    # ln[(n-m)!/(n+m)!] = lgamma(n-m+1) - lgamma(n+m+1)
    
    # Aplicando a exponencial em ambos os lados, sumimos com o ln na esquerda e retornamos ao que queríamos
    # (n-m)!/(n+m)! = exp(lgamma(n-m+1) - lgamma(n+m+1))
    ratio_factorials = exp(num - denom)

    # Trocando: "calcular dois números enormes e depois dividir" por "subtrair dois logs (que são números pequenos e seguros), depois aplicar exp"

    phase_factor = ComplexF64(im^(n - m)) # i^(n - m)
    amplitude_factor = (2.0 * n + 1.0)
    legendre_val = Plm(cos(α), n, m_abs) # Plm da biblioteca LegendrePolynomials

    return phase_factor * amplitude_factor * ratio_factorials * legendre_val
end

# Esta é a nova versão da sua função.
# Ela recebe a "tabela" de BSCs já pronta.
function partialwavexp(ψ₀, r, θ, ϕ, n_max, νi, bsc_precalculado, plm_precalculado)
    
    # Se estamos no centro
    if r == 0.0
        return νi == 0 ? ComplexF64(ψ₀) : ComplexF64(0.0)
    end

    k = 1000 * 2 * pi / 1.54
    kr = k * r
    ψ_total = ComplexF64(0.0)

    # O loop em 'm' desaparece pois m=νi.
    # O loop principal agora é só em 'n'.
    for n in abs(νi):n_max # n deve ser sempre >= |m|
        
        # Busca o BSC na tabela, em vez de recalcular!
        bsc = bsc_precalculado[n+1]

        legendre = plm_precalculado[n+1]
        
        # Pula o resto se o Plm for zero (ganho de performance)
        if legendre == 0.0
            continue
        end
        
        besselspher = sphericalbesselj(n, kr)
        ψ_total += bsc * besselspher * legendre
    end

    # A parte 'cis(m*ϕ)' da fórmula original agora fica fora do loop
    return ψ₀ * ψ_total * cis(νi * ϕ)
end

function makeplotpartial(nx, ny, rx, ry, ψ₀)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)

    k = 1000 * 2 * pi / 1.54

    n_max = 1000
    α = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ν = [0, 1, 5]

    mag = 10
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0

    heatmaps = Vector{Any}();

    for αi in α
        for νi in ν

            alfcount += 1 # Pois estava em 0 e o veotr começa em 1
            bessel_beam = Array{Float64}(undef, nx, ny) # Vetor que armazena as amplitudes do feixe de bessel para cada x,y

            # Essa parte vai pré-calcular os coeficiente de forma e os polinômios de legendre para cada valor de v, e alpha.
            bsc_precalculado = Vector{ComplexF64}(undef, n_max + 1)
            for n in 0:n_max
                bsc_precalculado[n+1] = bsccalc(n, νi, αi) 
            end

            plm_precalculado = Vector{Float64}(undef, n_max + 1)
            m_abs = abs(νi)
            for n in 0:n_max
                # Plm(0, n, m) é zero se n+m for ímpar.
                if n < m_abs || (n + m_abs) % 2 != 0
                    plm_precalculado[n+1] = 0.0
                else
                    plm_precalculado[n+1] = Plm(0.0, n, m_abs)
                end
            end

            @threads for j in 1:ny
                for i in 1:nx
                    xi = x[i]
                    yi = y[j]   

                    ρ = sqrt(xi^2 + yi^2)
                    r = sqrt(xi^2 + yi^2 + z^2)

                    if ρ == 0.0
                        ϕ = 0.0
                        # Não precisamos de θ, pois ele está "embutido" no plm_precalculado
                    else
                        ϕ = atan(yi, xi)
                    end
                    # if yi > 0
                    #     ϕ = acos(xi / ρ)
                    # else
                    #     if xi > 0
                    #         ϕ = 2 * pi - acos(xi / ρ)
                    #     else
                    #         ϕ = pi + acos(-xi / ρ)
                    #     end
                    # end 
                    θ = pi/2
                    #θ = acos(z/r);
                    # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                    bessel_beam[i, j] = abs(partialwavexp(ψ₀, r, θ, ϕ, n_max, νi, bsc_precalculado,plm_precalculado))^2
                end
            end
            push!(heatmaps, heatmap(1000 * x, 1000 * y, bessel_beam, xlabel=latexstring("x (mm)"), ylabel=latexstring("y (mm)"), title=latexstring('(' * alf[alfcount] * ')')))
            
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
    savefig(plot(heatmaps..., layout=(3, 3), size=(1200, 900)), "partial_wave_exp1.png")
end

# Mantenha suas funções bsccalc e partialwavexp como estão. Elas estão perfeitas.

# NOVA FUNÇÃO DE PLOT - Simples e correta
function makeplotpartial_displaced(nx, ny, rx, ry, ψ₀, x₀, y₀, z₀)
    
    # 1. Define o grid de visualização (exatamente como antes)
    z_grid = 0.0
    x_grid = range(-rx, rx, length=nx)
    y_grid = range(-ry, ry, length=ny)
    
    # 2. Define os parâmetros da simulação
    n_max = 1651 # Um valor razoável para um bom resultado 
    #n_max = Int64(ceil(kr + (big(405) / 100) * (kr^(1/3)) + 2))
    α = deg2rad(1.0) # Fixo em 1°, como na figura de referência
    ν_list = [0, 1, 5]

    alf = ['a', 'd', 'g']
    alfcount = 0
    heatmaps = Vector{Any}()

    # 3. Loop para cada subplot
    for νi in ν_list
        alfcount += 1
        bessel_beam = Array{Float64}(undef, nx, ny)

        # 4. PRÉ-CÁLCULO (A chave da velocidade)
        # Calculamos os BSCs e Plm para o feixe CENTRADO, uma única vez.
        # É rápido porque o somatório em 'm' não existe.
        bsc_precalculado = [bsccalc(n, νi, α) for n in 0:n_max]
        
        m_abs = abs(νi)
        plm_precalculado = [(n < m_abs || (n + m_abs) % 2 != 0) ? 0.0 : Plm(0.0, n, m_abs) for n in 0:n_max]
        
        println("Calculando para ν = $νi...")

        # 5. LOOP NOS PIXELS (com a lógica correta)
        @threads for j in 1:ny
            for i in 1:nx
                # Coordenadas do ponto do grid
                xi = x_grid[i]
                yi = y_grid[j]   

                # --- A LÓGICA SIMPLES E CORRETA ---
                # Em vez de usar (xi, yi) direto, calculamos as coordenadas
                # relativas ao centro do feixe (x₀, y₀)
                x_rel = xi - x₀
                y_rel = yi - y₀
                z_rel = z_grid - z₀ # A única linha nova!
                
                # Agora, convertemos essas coordenadas RELATIVAS para esféricas (r, θ, ϕ)
                r = sqrt(x_rel^2 + y_rel^2 + z_rel^2)
                ρ = sqrt(x_rel^2 + y_rel^2)
                
                ϕ = (ρ == 0.0) ? 0.0 : atan(y_rel, x_rel)
                θ = (r == 0.0) ? 0.0 : acos(z_rel/r) # Para z=0, θ é sempre pi/2
                
                # E chamamos a SUA função 'partialwavexp' com essas coordenadas relativas.
                bessel_beam[i, j] = abs2(partialwavexp(ψ₀, r, θ, ϕ, n_max, νi, bsc_precalculado, plm_precalculado))
            end
        end
        
        push!(heatmaps, heatmap(1000 * x_grid, 1000 * y_grid, bessel_beam, xlabel=latexstring("x_d (mm)"),ylabel=latexstring("y_d (mm)"), title=latexstring('(' * alf[alfcount] * ')'), margin=0Plots.mm))
    end

    savefig(plot(heatmaps..., layout=(3, 1), size=(500, 1300)), "feixe_deslocado_agora_sim.png")
    println("Simulação finalizada.")
end



# --- COMO USAR ---
# Deslocamento (x₀, y₀, z₀) = (100, 100, 1) em mm.
x_desloc = 0.1
y_desloc = 0.1
z_desloc = 0.01

# Pode rodar com uma resolução boa, porque este método é RÁPIDO.
@time makeplotpartial_displaced(200, 200, 0.2, 0.2, 1.0, x_desloc, y_desloc, z_desloc)
@time makeplotpartial_displaced(200, 200, 0.2, 0.2, 1.0, x_desloc, y_desloc, z_desloc*4)

#@time makeplotpartial(200, 200, 0.2, 0.2, 1)

