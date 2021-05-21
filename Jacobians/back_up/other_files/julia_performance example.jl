# Load Packages
using Optim
using Interpolations

# Define global constants
const carrying_capacity = 100
const r = 0.8
const discount_factor = 0.9
const convcrit = 1e-6
const ngrid = 1001
const stock_grid = range(0, carrying_capacity, length = ngrid)

# Define Functions
function stock_next(harvest, stock)
    new_stock = stock - harvest + r * stock .* (1 .- stock/carrying_capacity)

    if new_stock > 100.0
        new_stock = 100.0
    elseif new_stock < 0.0
        new_stock = 0.0
    end

    new_stock
end

# Something more complicated will evenutally go here
function profit(harvest)
    harvest
end

function sup_norm(a,b)
     maximum(abs.(a - b))
end

function check_convergence(a,b)
    sup_norm(a,b) < convcrit
end

function bellman(harvest, stock, valu_interp)
    current_reward = profit(harvest)
    future_stock = stock_next(harvest, stock)
    future_reward = valu_interp(future_stock)
    current_reward + discount_factor*future_reward
end

function valu_iteration()
    harvest = zeros(ngrid)
    valu = zeros(ngrid)
    valu_next = ones(ngrid)
    counter = 0
    converged = false

    while (converged == false) & (counter < 10000)

        harvest = zeros(ngrid)
        valu_next = zeros(ngrid)
        valu_interp = LinearInterpolation(stock_grid, valu)

        for idx in 1:ngrid
            stock = stock_grid[idx]
            objective = h -> - bellman(h, stock, valu_interp)
            result = optimize(objective, 0.0, stock)
            harvest[idx] = result.minimizer
            valu_next[idx] = -result.minimum
        end

        counter += 1  # Add 1 to the counter
        converged = check_convergence(valu,valu_next)
        valu = valu_next
    end
    harvest
end

valu_iteration() # Run once to compile
@time valu_iteration() # Time it.
