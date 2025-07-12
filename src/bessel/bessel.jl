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

# O algoritmo de miller usa recorrência reversa para evitar erros de arredondamento em "n" altos, pois no cálculo direto ordens altas erros se acumulam dando instabilidade numérica

# orderLim: até qual valor de n queremos calcularemos
# arg: valor de x em J_n(x)
# tol: até que ordem superior iremos fazer os cálculos para ter segurança e estabilidade
function besselj_miller(orderLim, arg, tol)
    # Assert: verifica se uma condição é verdadeira, se não for ele lança um erro com mensagem opcional
    @assert tol >= orderLim "Tolerance must be greater than the order to be calculated"

    N = tol
    values = Vector{typeof(arg)}(undef, N + 2) # armazena valores do tipo do argumento, não são inicializados com nada (undef) e tem tamanho N+2
    # values[1] = J0(x)

    # Chutes iniciais
    values[N+2] = 0 # O último índice é inicializado com 0 (arbitrário)
    values[N+1] = 1
    norm = 0.0

    for i in N+1:-1:2 # Começa em N+1, pois iremos precisar do seguinte e do anterior, se começasse em N+2, o seguinte extrapolaria os limites do vetor
        values[i-1] = (2 * (i - 1) / arg) * values[i] - values[i+1]
    end

    norm = values[1]
    for i in 3:2:N+2 # Apenas os i ímpares (ou seja, n pares) são multiplicados por 2
        norm += 2 * values[i];
    end

    norm = 1 / norm
    values = values[1:orderLim+1] # Truncamos nosso vetor, pois aí sim temos de 0 até ordemLim!
    return norm .* values # multiplicação elemento a elemento para normalizar todos os valores de values
end

# Função que calcula o feixe de bessel cilíndrico em coordenadas cilíndricas
function besbeam(ψ₀, α, k, order, ρ, ϕ, z) # ψ0, α, k, order = ν, ρ, ϕ, z      
    kz = k * cos(α) # Componente Longitudinal 
    kρ = k * sin(α) # Componente transversal
    return ψ₀ * besselj(order, kρ * ρ) * cis(order * ϕ) * cis(kz * z)
end # Feixes de bessel são descritos por funções de Bessel do primeiro tipo, são não difrativos e autocurativos

# A função a seguir, realiza a plotagem do gráfico
function makeplotbessel(nx, ny, rx, ry, ψ₀)

    z = 0 # Plot centrado em z = 0
    x = range(-rx, rx, nx) # Dentro dos limites (-rx,rx) criam-se "nx" pontos
    y = range(-ry, ry, ny)

    # Setando k
    # k = 1000 * 2 * big(pi) / big(1.54) # fixed wavelength 1.54mm (Maior precisão para k)
    k = 1000 * 2 * pi / 1.54 # TIREI OS BIGS (DESEMPENHO) AGORA FICOU (Float64)

    # Vamos simular para três tipos de ângulos áxicon α (1°, 10° e 40°)
    α = [deg2rad(1), deg2rad(10), deg2rad(40)]

    # As ordens dos feixes de Bessel que serão plotados
    ν = [0, 1, 5]

    # Fator de magnificação inicial para as coordenadas
    mag = 10

    # Nomeação dos subplots
    alf = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'] # 9 subplots
    alfcount = 0 # Contador para o array alf

    heatmaps = Vector{Any}()
    for αi in α
        for νi in ν

            alfcount += 1 # Pois estava em 0 e o veotr começa em 1
            bessel_beam = Array{Float64}(undef, nx, ny) # Vetor que armazena as amplitudes do feixe de bessel para cada x,y

            @threads for i in 1:nx
                # println("Thread ", Threads.threadid(), " está processando linha j=", j)
                for j in 1:ny
                    xi = x[i]
                    yi = y[j]   
                    ρ = sqrt(xi^2 + yi^2) # TIREI O BIG (DESEMPENHO)

                    if yi > 0
                        ϕ = acos(xi / ρ)
                    else
                        if xi > 0
                            ϕ = 2 * pi - acos(xi / ρ)
                        else
                            ϕ = pi + acos(-xi / ρ)
                        end
                    end

                    # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                    bessel_beam[j, i] = abs(besbeam(ψ₀, αi, k, νi, ρ, ϕ, z))^2
                end
            end
            push!(heatmaps, heatmap(1000 * x, 1000 * y, bessel_beam, xlabel=latexstring("x (mm)"), ylabel=latexstring("y (mm)"), title=latexstring('(' * alf[alfcount] * ')')))
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
    savefig(plot(heatmaps..., layout=(3, 3), size=(1200, 900)), "bessel_analytical.png")
end

function makeplotbessel_sideview(ny, nz, ry, rz, ψ₀)

    x = 0.0 # Plot centrado em x = 0

    # Dentro dos limites (-ry, ry) criam-se "ny" pontos e o mesmo para z
    y = range(-ry, ry, ny)
    z = range(-rz, rz, nz) 
    
    k = 1000 * 2 * pi / 1.54 # TIREI OS BIGS (DESEMPENHO) AGORA FICOU (Float64)

    # Vamos simular para três tipos de ângulos áxicon α (1°, 10° e 40°)
    α = [deg2rad(1), deg2rad(10), deg2rad(40)]

    # As ordens dos feixes de Bessel que serão plotados
    ν = [0, 1, 5]

    # Fator de magnificação inicial para as coordenadas
    mag = 10

    # Nomeação dos subplots
    alf = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'] # 9 subplots
    alfcount = 0 # Contador para o array alf

    heatmaps = Vector{Any}()
    for αi in α
        for νi in ν

            alfcount += 1 # Pois estava em 0 e o veotr começa em 1
            bessel_beam = Array{Float64}(undef, nz, ny) # Vetor que armazena as amplitudes do feixe de bessel para cada y,z

            @threads for i in 1:nz
                for j in 1:ny
                    zi = z[i]  
                    yi = y[j]
                    
                    # ρ = sqrt(xi^2 + yi^2)
                    ρ = abs(yi)

                    # Para evitar problemas de divisão por zero:
                    if ρ == 0.0 
                        ϕ = 0
                    end

                    # elseif yi > 0
                    #     ϕ = acos(xi / ρ)
                    # else
                    #     if xi >= 0 # Usar >= para incluir o caso xi=0, yi<0
                    #         ϕ = 2 * pi - acos(xi / ρ)
                    #     else
                    #         ϕ = pi + acos(-xi / ρ)
                    #     end
                    # end
                    # xi = x. Como x é fixado em zero. Podemos simplificar as próximas linhas em

                    ϕ = (yi >= 0) ? pi/2 : 3*pi/2 

                    # Agora vamos calcular a amplitude do feixe para cada z,y (ρ,ϕ)
                    # o Plots.jl espera que o array seja [length(y), length(z)] (vertical, horizontal).
                    bessel_beam[j, i] = abs(besbeam(ψ₀, αi, k, νi, ρ, ϕ, zi))^2
                end
            end
            push!(heatmaps, heatmap(1000 * z, 1000 * y, bessel_beam, xlabel=latexstring("z (mm)"), ylabel=latexstring("y (mm)"), title=latexstring('(' * alf[alfcount] * ')')))
        end
        y = y / mag
        z = z / mag
        mag = 4
    end
    savefig(plot(heatmaps..., layout=(3, 3), size=(1200, 900)), "bessel_sideview.png")
end

function makeplotbessel_desloc(nx, ny, rx, ry, ψ₀, x₀, y₀)

    z = 0 # Plot centrado em z = 0
    x = range(-rx, rx, nx) # Dentro dos limites (-rx,rx) criam-se "nx" pontos
    y = range(-ry, ry, ny)

    # Setando k
    # k = 1000 * 2 * big(pi) / big(1.54) # fixed wavelength 1.54mm (Maior precisão para k)
    k = 1000 * 2 * pi / 1.54 # TIREI OS BIGS (DESEMPENHO) AGORA FICOU (Float64)

    # Vamos simular para três tipos de ângulos áxicon α (1°, 10° e 40°)
    α = [deg2rad(1), deg2rad(10), deg2rad(40)]

    # As ordens dos feixes de Bessel que serão plotados
    ν = [0, 1, 5]

    # Fator de magnificação inicial para as coordenadas
    mag = 10

    # Nomeação dos subplots
    alf = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'] # 9 subplots
    alfcount = 0 # Contador para o array alf

    heatmaps = Vector{Any}()
    for αi in α
        for νi in ν

            alfcount += 1 # Pois estava em 0 e o veotr começa em 1
            bessel_beam = Array{Float64}(undef, nx, ny) # Vetor que armazena as amplitudes do feixe de bessel para cada x,y

            @threads for j in 1:ny
                yi = y[j]
                # println("Thread ", Threads.threadid(), " está processando linha j=", j)
                for i in 1:nx
                    xi = x[i]

                    x_rel = xi - x₀
                    y_rel = yi - y₀

                    ρ = sqrt(x_rel^2 + y_rel^2) 

                    ϕ = atan(y_rel, x_rel)

                    # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                    bessel_beam[i, j] = abs(besbeam(ψ₀, αi, k, νi, ρ, ϕ, z))^2
                end
            end
            push!(heatmaps, heatmap(1000 * x, 1000 * y, bessel_beam, xlabel=latexstring("x (mm)"), ylabel=latexstring("y (mm)"), title=latexstring('(' * alf[alfcount] * ')')))
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
    savefig(plot(heatmaps..., layout=(3, 3), size=(1200, 900)), "bessel_analytical_desloc.png")
end

@time makeplotbessel(200, 200, 0.2, 0.2, 1)
@time makeplotbessel_sideview(200, 200, 0.2, 0.2, 1)
@time makeplotbessel_desloc(200, 200, 0.2, 0.2, 1.0, -0.1, -0.1)

