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

function partialwavexp_deslocado(ψ₀, axiconang, order, k, r, θ, ϕ, x0, y0, z0, nmax)
    kr = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = atan(y0, x0)
    rho0 = sqrt(x0^2 + y0^2)
    
    psi = 0.0 + 0.0im

    # Pré-calcula os fatores de deslocamento para cada m
    for n in 0:nmax
        spher = sphericalbesselj(n, kr)
        for m in -n:n
            # Fatores de deslocamento
            bessel_factor = besselj(m - order, krho * rho0)
            phase_factor = cis(-(m - order) * phi0) * cis(-kz * z0)
            # BSC centrado
            BSC = bsccalc(n, m, axiconang)
            # Soma todos os modos m
            psi += phase_factor * bessel_factor * BSC * spher * Plm(cos(θ), n, m) * cis(m * ϕ)
        end
    end
    return ψ₀ * psi
end

function makeplotpartial(nx, ny, rx, ry, ψ₀, x0, y0, z0)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)

    k = 1000 * 2 * pi / 1.54

    n_max = 100
    # nmax = Int(ceil(kr + (405 / 100) * (kr^(1/3)) + 2))
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
                   
                    θ = pi/2
                    #θ = acos(z/r);
                    # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                    bessel_beam[i, j] = abs(partialwavexp_deslocado(ψ₀, αi, νi, k, r, θ, ϕ, x0, y0, z0,n_max))^2
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

@time makeplotpartial(200, 200, 0.2, 0.2, 1, 0.0, 0.0, 0.0)

