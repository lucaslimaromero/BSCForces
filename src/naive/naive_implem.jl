using Bessels
using SpecialFunctions
using LegendrePolynomials
using Plots
using BenchmarkTools

global fact = zeros(typeof(big(1)*big(1)), 10000)
fact[1] = 1

function besselN(N, z)
    cosz = cos(z)
    sinz = sin(z)
    frac = (4 * N^2 - 1) / (8 * z)
    nfrac = (4 * (N+1)^2 - 1) / (8 * z)
    if N % 2 == 0 
        return (-sinz + cosz +  frac * (cosz + sinz)) / (-cosz - sinz + nfrac * (-sinz + cosz))
    else
        return (sinz + cosz +  frac * (-cosz + sinz)) / (cosz - sinz + nfrac * (sinz + cosz))
    end
end

function besselj(orderLim, arg, tol)
    N = tol * orderLim
    values = Vector{typeof(arg)}(undef, N+2)    
    values[N+1] = besselN(N, arg)
    values[N+2] = 1
    norm = 0.0
    for i in N+1:-1:2
        values[i-1] = (2 * (i-1) / arg) * values[i] - values[i+1]
    end 
    norm = values[1]
    for i in 3:2:N+2
        norm += 2 * values[i]
    end
    norm = 1/norm
    values = values[1:orderLim+1]
    return norm.*values
end

function factorial(x::T) where T<:BigInt
    if x == 0
        return 1
    end
    if fact[x] == zero(BigInt)
        i::typeof(x*x) = x
        while fact[i] == zero(typeof(x*x))
            i-=1
        end
        i+=1
        while i != x + 1
            fact[i] = fact[i-1] * i
            i+=1
        end
    end
    return fact[x]
end

function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = 1000 * 2 * big(pi) / big(1.54) #fixed wavelength 1.54mm
    kz = k * big(cos(axiconang))
    krho = k * big(sin(axiconang))
    return big(psiamp * SpecialFunctions.besselj(order, krho * rho) * cis(order * phi) * cis(kz * z))
end

function bsccalc(n, m, axiconang, order, krho, kz, phi0, z0, rho0)
    fract = (im ^ (n - m)) * (2 * n + 1) * factorial(big(n - m)) / factorial(big(n+m))
    special = Plm(cos(axiconang), n, m)
    return fract * special 
end

function partialwavexp(psiamp, axiconang, order, r, theta, phi)
    k = 1000 * 2 * big(pi) / big(1.54)
    kr = k * r
    norm = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = z0 = rho0 = 0
    nmax = 30 #Int64(ceil(norm + (big(405) / 100) * (norm^(1/3)) + 2))
    psi = 0
    besselm = []
    cisvalm = []
    for m in -nmax:nmax
        push!(cisvalm,  cis(-(m + order) * phi0) * cis(-kz * z0))
        push!(besselm, SpecialFunctions.besselj(m - order, krho * rho0))
    end
    for n in 0:nmax
        spher = SpecialFunctions.sphericalbesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, axiconang, order, krho, kz, phi0, z0, rho0)
            psi += cisvalm[m + nmax + 1] * besselm[m + nmax + 1] * BSC * spher * Plm(cos(theta), n, m) * cis(m * phi)
        end
    end
    return psi * psiamp
end


function makeplotbes(nx, ny, ry, rx, psiamp)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)
    axang = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ord = [0, 1, 5]
    mag = 10
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0
    htmaps = Vector{Any}()
    for axiconang in axang
        for order in ord
            beamaux = Vector{BigFloat}()
            alfcount += 1
            for j in y, i in x
                rho =  sqrt(big(i)^2 + big(j)^2)
                if j > 0
                    phi = big(acos(big(i) / rho))
                else
                    if i > 0
                        phi = 2 * pi - big(acos(big(i) / rho))
                    else
                        phi = pi +  big(acos(-big(i) / rho))
                    end
                end
                push!(beamaux, abs(besbeam(psiamp, big(axiconang), order, rho, phi, z)))
            end
            push!(htmaps, heatmap(x, y,  reshape(beamaux, (nx, ny)), xlabel ='x', ylabel='y', title = '(' * alf[alfcount] * ')' , titlefontsize = 12, xlabelfontsize = 9, ylabelfontsize = 9))
            println(typeof(reshape(beamaux, (nx, ny))))
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
    savefig(plot(htmaps..., layout = (3, 3),  xtickfontsize = 9, yfontsize = 9, size =(1200,900)), "plotorg.png")
end



function makepltpartial(nx, ny, ry, rx, psiamp)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)
    axang = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ord = [0, 1, 5]
    mag = 10
    testc = 0
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0
    htmaps = Vector{Any}()
    for axiconang in axang
        for order in ord
            beamaux = Vector{BigFloat}()
            alfcount+=1
            for j in y, i in x
                rho =  sqrt(big(i)^2 + big(j)^2)
                r = sqrt(big(i)^2 + big(j)^2 + z^2)
                if j > 0
                    phi = big(acos(big(i) / rho))
                else
                    if i > 0
                        phi = 2 * pi - big(acos(big(i) / rho))
                    else
                        phi = pi +  big(acos(-big(i) / rho))
                    end
                end
                theta = acos(z)
                push!(beamaux,abs(partialwavexp(psiamp, big(axiconang), order, r, theta, phi)))
                println(testc)
                testc+=1
            end
            push!(htmaps, heatmap(x, y,  reshape(beamaux, (nx, ny)), xlabel ='x', ylabel='y', title = '(' * alf[alfcount] * ')' , titlefontsize = 12, xlabelfontsize = 9, ylabelfontsize = 9))
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
savefig(plot(htmaps..., layout = (3, 3),  xtickfontsize = 9, yfontsize = 9, size =(1200,900)), "plotpart_test.png")
end

makepltpartial(200, 200, 0.2, 0.2, 1)
# makeplotbes(200, 200, 0.2, 0.2, 1)







