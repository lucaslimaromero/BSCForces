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

function partialwavexplight(psiamp, axiconang, order, r, θ, ϕ)
    k = 1000 * 2 * pi / 1.54
    kr = k * r
    
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    ϕ0 = z0 = rho0 = 0
    nmax = 150 #Int64(ceil(kr + (big(405) / 100) * (kr^(1/3)) + 2))
    # psi = 0

    psi = ComplexF64(0.0)

    besselm = Vector{Float64}(undef, 2 * nmax + 1)
    cisvalm = Vector{ComplexF64}(undef, 2 * nmax + 1)

    for m in -nmax:nmax
        idx = m + nmax + 1
        cisvalm[idx] = cis(-(m + order) * ϕ0) * cis(-kz * z0)
        besselm[idx] = besselj(m - order, krho * rho0)
    end

    for n in 0:nmax
        spher = SpecialFunctions.sphericalbesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, axiconang)
            psi += cisvalm[m + nmax + 1] * besselm[m + nmax + 1] * BSC * spher * Plm(cos(θ), n, m) * cis(m * ϕ)
        end
    end
    return psi * psiamp
end

function makeplotpartial(nx, ny, rx, ry, ψ₀)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)

    k = 1000 * 2 * pi / 1.54

    α = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ν = [0, 1, 5]
    mag = 10
    testc = 0
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0

    heatmaps = Vector{Any}(undef, 9);
    for αi in α
        for νi in ν

            alfcount += 1 # Pois estava em 0 e o veotr começa em 1
            bessel_beam = Array{Float64}(undef, nx, ny) # Vetor que armazena as amplitudes do feixe de bessel para cada x,y

            @threads for j in 1:ny
                yi = y[j]
                for i in 1:nx
                    xi = x[i]

                    ρ = sqrt(xi^2 + yi^2)
                    r = sqrt(xi^2 + yi^2 + z^2)

                    if yi > 0
                        ϕ = acos(xi / ρ)
                    else
                        if xi > 0
                            ϕ = 2 * pi - acos(xi / ρ)
                        else
                            ϕ = pi + acos(-xi / ρ)
                        end
                    end # z = cosθ

                    θ = acos(z);
                    # Agora vamos calcular a amplitude do feixe para cada x,y (ρ,ϕ)
                    bessel_beam[i, j] = abs(partialwavexplight(ψ₀, αi, νi, r, θ, ϕ))^2
                    # println("Thread ", Threads.threadid(), " i=", i, " j=", j, " xi=", xi, " yi=", yi)
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

@time makeplotpartial(30, 30, 0.2, 0.2, 1)

