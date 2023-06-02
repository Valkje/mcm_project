println("Importing packages...")

using LinearAlgebra
using Distributions
using Plots

# Make image directory
im_path = "images"
mkpath(im_path)

######### (a) #########

println("Running (a)...")

function quadraticpotential_grad(q)
    q
end

function leapfrog(q0, p0, tau, l, mass, potential_grad::Function)
    # Determine dimensionality
    d = length(q0)

    # Also accomodate initial values
    qs = [Vector{Float64}(undef, d) for _ = 1:(l+1)]
    ps = [Vector{Float64}(undef, d) for _ = 1:(l+1)]

    qs[1] = q0
    ps[1] = p0

    for i = 2:(l+1)
        p_half = ps[i - 1] .- tau / 2 .* potential_grad(qs[i - 1])
        qs[i] = qs[i - 1] .+ tau * inv(mass) * p_half
        ps[i] = p_half .- tau / 2 .* potential_grad(qs[i])
    end

    # Discard initial value
    (qs[2:(l+1)], ps[2:(l+1)])
end

######### (b) #########

println("Running (b)...")

q0 = [0.]
p0 = [1.]
tau = 1.2
l = 20
mass = [1.;;]

(ps1, qs1) = leapfrog(q0, p0, tau, l, mass, quadraticpotential_grad)

tau = 0.3
l = 80

(ps2, qs2) = leapfrog(q0, p0, tau, l, mass, quadraticpotential_grad)

qs1_flat = [v[1] for v in qs1]
ps1_flat = [v[1] for v in ps1]

qs2_flat = [v[1] for v in qs2]
ps2_flat = [v[1] for v in ps2]

t = collect(0:0.01:2π)
q_exact = -cos.(t .+ π/2)
p_exact = sin.(t .+ π/2)

plot(q_exact, p_exact)
plot!(qs1_flat, ps1_flat)

savefig(joinpath(im_path, "leapfrog_large_tau.png"))

t = collect(0:0.01:2π)
q_exact = -cos.(t .+ π/2)
p_exact = sin.(t .+ π/2)

plot(q_exact, p_exact)
plot!(qs2_flat, ps2_flat)

savefig(joinpath(im_path, "leapfrog_small_tau.png"))

######### (f) #########

println("Running (f)...")

function quadraticpotential(q)
    0.5 * q' * q
end

function hamiltonian(q, p, mass, potential::Function)
    potential(q) + 0.5 * p' * inv(mass) * p
end

function hmc(n, q0, tau_dist, l, potential::Function, potential_grad::Function)
    # Determine dimensionality
    d = length(q0)

    # Also accomodate initial values
    qs = [Vector{Float64}(undef, d) for _ = 1:(n+1)]
    ps = [Vector{Float64}(undef, d) for _ = 1:(n+1)]

    # Julia can infer the size of I from the mathematical context
    mass = I
    norm = MvNormal(zeros(Float64, d), mass)

    qs[1] = q0

    for t = 2:(n+1)
        ps[t - 1] = vec(rand(norm, 1))

        tau = rand(tau_dist, 1)[1]
        # Note that qs_star has a different dimensionality (l) than qs (n)
        (qs_star, ps_star) = leapfrog(qs[t - 1], ps[t - 1], tau, l, mass, potential_grad)
        q_star = qs_star[l]
        p_star = ps_star[l]

        h_prev = hamiltonian(qs[t - 1], ps[t - 1], mass, potential)
        h_star = hamiltonian(q_star, p_star, mass, potential)

        alpha = min(1, exp(h_prev - h_star))

        if rand() < alpha
            qs[t] = q_star
            ps[t] = p_star
        else
            qs[t] = qs[t - 1]
            ps[t] = ps[t - 1]
        end
    end

    (qs[2:(n+1)], ps[2:(n+1)])
end

tau_dist = Uniform(0.5, 1.5)
samples, _ = hmc(100000, q0, tau_dist, l, quadraticpotential, quadraticpotential_grad)

samples_flat = [s[1] for s in samples]
histogram(samples_flat)

x = -5:0.1:5
plot!(x, pdf.(Normal(0, 1), x) .* 1e4, linewidth = 3)

savefig(joinpath(im_path, "hmc_test.png"))

######### (g) #########

println("Running (g)...")

function rw_mh(n, start, prop_dist, target::Function)
    d = length(start)

    samples = [Vector{Float64}(undef, d) for _ = 1:(n+1)]
    samples[1] = start

    for i = 2:(n+1)
        sample_prop = vec(rand(prop_dist, 1)) + samples[i - 1]

        dens_prev = target(samples[i - 1])
        dens_prop = target(sample_prop)
        alpha = min(1, dens_prop / dens_prev)

        if rand() < alpha
            samples[i] = sample_prop
        else
            samples[i] = samples[i - 1]
        end
    end

    samples[2:(n+1)]
end

target(x) = exp(-quadraticpotential(x))
samples = rw_mh(1000000, [0.], MvNormal(zeros(1), I), target)
histogram([s[1] for s in samples])

x = -5:0.1:5
plot!(x, target.(x) .* 2e4, linewidth = 3)

savefig(joinpath(im_path, "mh_test.png"))

######### (h) #########

println("Running (h)...")

function rand_sigma_mh(n, start, sigma_dist, target::Function)
    mh_samples = [Vector{Float64}(undef, d) for _ = 1:(n+1)]
    mh_samples[1] = start

    mu = zeros(size(start))

    for i = 2:(n+1)
        sigma = rand(sigma_dist, 1)[1]
        prop_dist = MvNormal(mu, sigma * I)

        mh_samples[i] = rw_mh(1, mh_samples[i - 1], prop_dist, target)[1]
    end

    mh_samples[2:(n+1)]
end

sigmas = collect(0.01:0.01:1)
Σ = Diagonal(sigmas)
d = size(Σ, 1)
μ = zeros(d)
target(x) = pdf(MvNormal(μ, Σ), x)

n = 1000
l = 150
tau_dist = Uniform(0.0104, 0.0156)
potential(x) = -log(target(x))
potential_grad(x) = @. (1 / sigmas) * x
samples_hmc, _ = hmc(n, μ, tau_dist, l, potential, potential_grad)

n = 150000
sigma_dist = Uniform(0.0176, 0.0264)
samples_mh = rand_sigma_mh(n, μ, sigma_dist, target)
samples_mh_thinned = samples_mh[1:150:length(samples_mh)];

# Plotting
last_comp_hmc = [v[d] for v in samples_hmc]
last_comp_mh_th = [v[d] for v in samples_mh_thinned]

plot(last_comp_hmc, labels = "HMC")
p1 = plot!(last_comp_mh_th, labels = "MH (thinned)")

p2 = plot(last_comp_hmc, legend = false)

p3 = plot(last_comp_mh_th, color = 2, legend = false)

plot(p1, p2, p3, layout = (3, 1))

savefig(joinpath(im_path, "sample_time_course.png"))

# Turns the vector of vectors inside-out (components outside, samples inside)
comps_hmc = [[v[i] for v in samples_hmc] for i = 1:d]
means_hmc = mean.(comps_hmc)
plot(sigmas, means_hmc, labels = "HMC")

comps_mh = [[v[i] for v in samples_mh] for i = 1:d]
means_mh = mean.(comps_mh)
p1 = plot!(sigmas, means_mh, labels = "MH (thinned)")

p2 = plot(sigmas, means_hmc, legend = false)

p3 = plot(sigmas, means_mh, color = 2, legend = false)

plot(p1, p2, p3, layout = (3, 1))

savefig(joinpath(im_path, "sample_mean_vs_sd.png"))

println("Done!")
