#using Bessels # Alta precisão para funções de Bessel
using SpecialFunctions # inclui funções como besselj e sphericalbesselj
using LegendrePolynomials
using Plots
using BenchmarkTools
using Test 
using LaTeXStrings
using Plots.PlotMeasures

# α: axiconang
# h = kρ = k * sinα      
# β = kz = k * cosα      

# Armazenar fatoriais previamente calculados para evitar recálculos
global fact = zeros(typeof(big(1)*big(1)), 10000) # Preenchemos um vetor global de 10000 posições com zeros
fact[1] = 1 # Ex.: fact[5] deverá armazenar 120

function factorial(x::T) where T<:BigInt
    if x == 0
        return 1 # 0! = 1, Caso Base
    end

    if fact[x] == zero(BigInt) # Se o fatorial de x não foi calculado, calcularemos recursivamente
        i::typeof(x*x) = x # i = x
        while fact[i] == zero(typeof(x*x)) # Procura o maior fatorial i já calculado antes de x
            i -= 1
        end
        # i = 5, p ex e queremos achar x! = 7! daí decrescemos i = 6, vimos que fact[6] era zero, não tinha sido calculados e decrescemos novamente
        i += 1 # Vamos para i = 6, pois o de 5 já temos (step back) e calculamos 6!
        while i != x + 1 # (x + 1 = 8)
            fact[i] = fact[i-1] * i # 6! = 6.5!
            i += 1 # i = 7
        end
    end

    return fact[x]
end

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

# Função que calcula os Beam Shape Coefficients (BSCs)
function bsccalc(n, m, α) 
    fract = (im^(n - m)) * (2 * n + 1) * fact(big(n - m)) / fact(big(n + m))
    legendre = Plm(cos(α), n, m)
    return fract * legendre
end # Esssa função só calcula a primeira parcela do BSC, ignora a J, exp de ϕ0 e exp de z0 (calculado e multiplicado depois)
# As simulações são feitas para z0 = 0, portanto a dependência exponencial de z0 é 1.

# Função que calcula o feixe de bessel cilíndrico em coordenadas cilíndricas
function besbeam(ψ₀, α, k, order, ρ, ϕ, z) # ψ0, α, k, order = ν, ρ, ϕ, z      
    kz = k * cos(α) # Componente Longitudinal 
    kρ = k * sin(α) # Componente transversal
    return ψ₀ * besselj(order, kρ * ρ) * cis(order * ϕ) * cis(kz * z)
end # Feixes de bessel são descritos por funções de Bessel do primeiro tipo, são não difrativos e autocurativos

# Calcula o campo acústico completo (3D) simulando a propagação de um feixe de Bessel em todo o volume definido por x,y,z
function completebeam(x, y, z, α, order, amp, k, x0, y0, z0)

    # Array vazio para armazenar o campo
    beam = []

    for zi in z, yi in y, xi in x
        # Calcularemos a contribuição em cada ponto usando partialwavexp
        push!(beam, partialwavexp(amp, α, order, k, 0, 0, 0, xi - x0, yi - y0, zi - z0))
    end

    return reshape(beam, (size(x)[1], size(y)[1], size(z)[1]))
end 

# Essa função simplesmente calcula a amplitude da onda parcial dado um ponto no espaço
function partialwavexp(ψ₀, α, order, k, r, θ, ϕ, x0, y0, z0)

    kr = k * r
    kρ = k * sin(α)
    kz = k * cos(α)

    # Para achar os fatores de forma generalizados irei precisar de ρ0, ϕ0 e z0.
    # Consigo obter esses valores a partir da relação entre as coordenadas cilíndricas e as cartesianas x0, y0, z0 fornecidas 
    ϕ0 = Base.atan(y0/x0)
    ρ0 = (x0^2 + y0^2)^(1/2)
    # z0 = z0

    # Agora, tenhamos uma boa perspectiva de nmax, já que o computador não resolve somatórios infinitos, precisaremos truncar inteligentemente
    nmax = Int64(ceil(kr + 4.05 * (kr)^(1/3) + 2))

    ψ = 0
    besselm = []
    expm = []

    # Equação (9) do meu projeto, veja que calculamos a parte que depende de m dos BSCs
    for m in -nmax:nmax
        push!(expm, cis(-(m - order) * ϕ0) * cis(-kz * z0))
        push!(besselm, besselj(m - order, Float64(kρ * ρ0)))
    end

    for n in 0:nmax
        besselspher = sphericalbesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, α) 
            ψ += BSC * besselm[n + nmax + 1] * expm[m + nmax + 1] * cis(m * ϕ) * besselspher * Plm(cos(θ), n, m)
        end # Só calculamos legendre para os BSCs, e não contemplamos o legendre de fora
    end

    return ψ * ψ0
end

# println(besselj0(3.0))
# println(besselj_miller(0,3.0,2))
# println("Pressione Enter para fechar...")
# readline()

# A função a seguir, realiza a plotagem do gráfico
function makeplotbessel(nx, ny, rx, ry, ψ₀)
    z = 0 # Plot centrado em z = 0
    x = range(-rx, rx, nx) # Dentro dos limites (-rx,rx) criam-se "nx" pontos
    y = range(-ry, ry, ny)

    # Setando k
    k = 1000 * 2 * big(pi) / big(1.54) # fixed wavelength 1.54mm

    # Vamos simular para três tipos de ângulos áxicon α
    α = [deg2rad(1), deg2rad(10), deg2rad(40)]

    # As ordens dos feixes de Bessel que serão plotados
    ν = [0, 1, 5]

    # Fator de magnificação inicial para as coordenadas
    mag = 10

    # Nomeação dos subplots
    alf = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ] # 9 subplots
    alfcount = 0 # Contador para o array alf

    heatmaps = Vector{Any}()
     
    for αi in α # Cada linha tem um ângulo de áxicon
        for νi in ν # Cada coluna é uma ordem

            bessel_beam = Vector{BigFloat}() # Vetor que armazena as amplitudes do feixe de bessel para cada x,y
            # Depois ele será convertido em uma matriz (nx, ny)
            alfcount += 1 # Pois estava em 0 e o veotr começa em 1

            for yi in y, xi in x

                ρ = sqrt(big(xi)^2 + big(yi)^2) 

                if yi > 0 # Se está no 1° ou segundo quadrante (yi > 0)
                    ϕ = big(acos(big(xi) / ρ))
                else
                    if xi > 0
                        ϕ = 2 * pi - big(acos(big(xi) / ρ))
                    else
                        ϕ = pi + big(acos(- big(xi) / ρ))
                    end
                end

                # # Testar
                # ϕ = atan(big(yi), big(xi))
                # if ϕ < 0
                #     ϕ += 2 * pi
                # end

                # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                push!(bessel_beam, abs(besbeam(ψ₀, αi, k, νi, ρ, ϕ, z))^2)
            end
            push!(heatmaps, heatmap(1000 * x, 1000 * y, reshape(bessel_beam, (nx, ny)), xlabel= latexstring("x (mm)"), ylabel= latexstring("y (mm)"), title = latexstring('(' * alf[alfcount] * ')'), titlefontsize = 16, xlabelfontsize = 12, ylabelfontsize = 12, xtickfontsize=9, ytickfontsize=9, left_margin=5mm, bottom_margin = 5mm))
        end

        # Após cada linha de subplots
        x = x / mag
        y = y / mag 
        mag = 4
    end
    savefig(plot(heatmaps..., layout = (3, 3), xtickfontsize = 9, yfontsize = 9, size = (1200,900)), "plot.png")
end


@time makeplotbessel(200, 200, 0.2, 0.2, 1)