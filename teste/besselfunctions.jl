using SpecialFunctions
using Plots
using BenchmarkTools
include("../src/naive/mynaive.jl")

# Configuração dos parâmetros
nx = 200
ny = 200
range_lim = 1.0
order = 1
miller_tolerance = 22

# Criando grids
x_vals = range(0, range_lim, length=nx)
y_vals = range(0, range_lim, length=ny)

# Pré-alocação das matrizes com tipo explícito
j1_matrix = zeros(Float64, ny, nx)
j1_miller = zeros(Float64, ny, nx)

# Plot dos erros
j1_error = zeros(Float64, ny, nx)

# Preenchendo as matrizes
for j in 1:ny
    y = y_vals[j]
    for i in 1:nx
        x = x_vals[i]
        rho = sqrt(x^2 + y^2)
        
        j1_matrix[j, i] = besselj(order, rho)  # Usando besselj em vez de besselj1 para consistência
        
        # Convertendo para Float64 diretamente
        miller_results = besselj_miller(order, rho, miller_tolerance)
        j1_miller[j, i] = Float64(miller_results[order + 1])
        

        error = abs(j1_miller[j, i] - j1_matrix[j, i]) / abs(j1_matrix[j, i])
        min_error = 1e-15  # Valor mínimo para evitar log10(0)
        relative_error = max(100 * error, min_error)
        j1_error[j, i] = log10(relative_error)
    end
end

# Criando os heatmaps
p1 = heatmap(x_vals, y_vals, j1_matrix,
            aspect_ratio=1,
            xlabel="x",
            ylabel="y",
            title="J₁(ρ) Nativo",
            titlefontsize=12,
            xlims=(0, range_lim), ylims=(0, range_lim))

p2 = heatmap(x_vals, y_vals, j1_miller,
            aspect_ratio=1,
            xlabel="x",
            ylabel="y",
            title="J₁(ρ) por Miller (N=22)",
            titlefontsize=12,
            xlims=(0, range_lim), ylims=(0, range_lim))

p3 = heatmap(x_vals, y_vals, j1_error,
            aspect_ratio=1,
            xlabel="x",
            ylabel="y",
            title="Erro Relativo",
            titlefontsize=12,
            xlims=(0, range_lim), ylims=(0, range_lim),
            clims=(-15, -11))

# Layout combinado com 3 colunas
plot_layout = @layout [a b c]
p_combined = plot(p1, p2, p3, layout=plot_layout, size=(1200, 400))

# Salvando e exibindo
savefig(p_combined, "bessel_comparison_subplot.png")
display(p_combined)

println("Pressione Enter para fechar...")
readline()