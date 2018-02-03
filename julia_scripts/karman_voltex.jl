using Plots
using DataFrames
using CSV

# Constants and Parameters
MX = 401
I1 = 95
I2 = 105
MY = 201
J1 = 95
J2 = 105
MAXITP = 10000
OMEGAP = 1.0
ERRORP = 1e-4
MAXSTEP = 5000
IO_FREQUENCY = 10
Re = 70

# Functions
function BoundaryConditionforP(p::Array{Float64, 2})
    for j in 1:MY
        p[1, j] = 0.0
        p[MX, j] = 0.0
    end

    for i in 1:MX
        p[i, 1] = 0.0
        p[i, MY] = 0.0
    end

    p[I1, J1] = p[I1-1, J1-1]
    p[I1, J2] = p[I1-1, J2+1]
    p[I2, J1] = p[I2+1, J1-1]
    p[I2, J2] = p[I2+1, J2+1]

    for j in J1+1:J2-1
        p[I1, j] = p[I1-1, j]
        p[I2, j] = p[I2+1, j]
    end

    for i in I1+1:I2-1
        p[i, J1] = p[i, J1-1]
        p[i, J2] = p[i, J2+1]
    end
end

function BoundaryConditionforV(u::Array{Float64, 2}, v::Array{Float64, 2})
    for j in 1:MY
        u[1, j] = 1.0
        v[1, j] = 0.0
        u[MX, j] = 2.0 * u[MX-1, j] - u[MX-2, j]
        v[MX, j] = 2.0 * v[MX-1, j] - v[MX-2, j]
    end

    for i in 1:MX
        u[i, 1] = 2.0 * u[i, 2] - u[i, 3]
        v[i, 1] = 2.0 * v[i, 2] - v[i, 3]
        u[i, MY] = 2.0 * u[i, MY-1] - u[i, MY-2]
        v[i, MY] = 2.0 * v[i, MY-1] - v[i, MY-2]
    end

    for i in I1:I2
        for j in J1:J2
            u[i, j] = 0.0
            v[i, j] = 0.0
        end
    end
end

function PoissonEquation(u::Array{Float64, 2}, v::Array{Float64, 2}, p::Array{Float64, 2},
                         dx::Float64, dy::Float64, dt::Float64)
    rhs = zeros(MX, MY)

    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 <= i <= I2) && (J1 <= j <= J2)
            else
                ux = (u[i+1, j] - u[i-1, j]) / (2.0*dx)
                uy = (u[i, j+1] - u[i, j-1]) / (2.0*dy)
                vx = (v[i+1, j] - v[i-1, j]) / (2.0*dx)
                vy = (v[i, j+1] - v[i, j-1]) / (2.0*dy)
                rhs[i, j] = (ux + vy) / dt - (ux^2 + 2.0 * uy * vx + vy^2)
            end
        end
    end

    @inbounds for iteration in 1:MAXITP
        res = 0.0
        for i in 2:MX-1
            for j in 2:MY-1
                if (I1 <= i <= I2) && (J1 <= j <= J2)
                else
                    dp = (p[i+1, j] + p[i-1, j]) / (dx^2) + (p[i, j+1] + p[i, j-1]) / (dy^2) -rhs[i, j]
                    dp = dp / (2.0 / (dx^2) + 2.0 /(dy^2)) - p[i, j]
                    res += dp^2
                    p[i, j] += OMEGAP * dp
                end
            end
        end
        BoundaryConditionforP(p)
        res = sqrt(res / (MX * MY))
        if res < ERRORP
            resp = res
            itrp = iteration
            break
        end
    end
end

function VelocityEquation(u::Array{Float64, 2}, v::Array{Float64, 2}, p::Array{Float64, 2},
                          dx::Float64, dy::Float64, dt::Float64)
    urhs = zeros(MX, MY)
    vrhs = zeros(MX, MY)
    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 <= i <= I2) && (J1 <= j <= J2)
            else
                urhs[i, j] = -(p[i+1, j] - p[i-1, j]) / (2.0*dx)
                vrhs[i, j] = -(p[i, j+1] - p[i, j-1]) / (2.0*dy)
            end
        end
    end

    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 <= i <= I2) && (J1 <= j <= J2)
            else
                urhs[i, j] += (u[i+1, j] - 2.0*u[i, j] + u[i-1, j]) / (Re*dx^2) + (u[i, j+1] - 2.0*u[i, j] + u[i, j-1]) / (Re*dy^2)
                vrhs[i, j] += (v[i+1, j] - 2.0*v[i, j] + v[i-1, j]) / (Re*dx^2) + (v[i, j+1] - 2.0*v[i, j] + v[i, j-1]) / (Re*dy^2)
            end
        end
    end

    for j in J1+1:J2-1
        u[I1+1, j] = 2.0*u[I1, j] - u[I1-1, j]
        u[I2-1, j] = 2.0*u[I2, j] - u[I2+1, j]
        v[I1+1, j] = 2.0*v[I1, j] - v[I1-1, j]
        v[I2-1, j] = 2.0*v[I2, j] - v[I2+1, j]
    end

    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 < i < I2) && (J1 < j < J2)
            elseif i == 2
                urhs[i, j] -= u[i, j] * (-u[i+2, j] + 8.0*(u[i+1, j] - u[i-1, j]) + u[i-1, j]) / (12.0*dx) +
                abs(u[i, j]) * (u[i+2, j] - 4.0*u[i+1, j] + 6.0*u[i, j] - 4.0*u[i-1, j] + u[i-1, j]) / (4.0*dx)
                vrhs[i, j] -= u[i, j] * (-v[i+2, j] + 8.0*(v[i+1, j] - v[i-1, j]) + v[i-1, j]) / (12.0*dx) +
                abs(u[i, j]) * (v[i+2, j] - 4.0*v[i+1, j] + 6.0*v[i, j] - 4.0*v[i-1, j] + v[i-1, j]) / (4.0*dx)
            elseif i == MX-1
                urhs[i, j] -= u[i, j] * (-2.0 * u[i+1, j] + u[i, j] + 8 * (u[i+1, j] - u[i-1, j]) + u[i-2, j]) / (12.0*dx) +
                abs(u[i, j]) * (2.0 * u[i+1, j] - u[i, j] - 4.0 * u[i+1, j] + 6.0 * u[i, j] - 4.0 * u[i-1, j] + u[i-2, j]) / (4.0*dx)
                vrhs[i, j] -= u[i, j] * (-2.0 * v[i+1, j] + v[i, j] + 8 * (v[i+1, j] - v[i-1, j]) + v[i-2, j]) / (12.0*dx) +
                abs(u[i, j]) * (2.0 * v[i+1, j] - v[i, j] - 4.0 * v[i+1, j] + 6.0 * v[i, j] - 4.0 * v[i-1, j] + v[i-2, j]) / (4.0*dx)
            else
                urhs[i, j] -= u[i, j] * (-u[i+2, j] + 8.0 * (u[i+1, j] - u[i-1, j]) + u[i-2, j]) / (12.0*dx) +
                abs(u[i, j]) * (u[i+2, j] - 4.0 * u[i+1, j] + 6.0 * u[i, j] - 4.0 * u[i-1, j] + u[i-2, j]) / (4.0*dx)
                vrhs[i, j] -= u[i, j] * (-v[i+2, j] + 8.0 * (v[i+1, j] - v[i-1, j]) + v[i-2, j]) / (12.0*dx) +
                abs(u[i, j]) * (v[i+2, j] - 4.0 * v[i+1, j] + 6.0 * v[i, j] - 4.0 * v[i-1, j] + v[i-2, j]) / (4.0*dx)
            end
        end
    end

    for i in I1+1:I2-1
        u[i, J1 + 1] = 2.0 * u[i, J1] - u[i, J1 - 1]
        u[i, J2 - 1] = 2.0 * u[i, J2] - u[i, J2 + 1]
        v[i, J1 + 1] = 2.0 * v[i, J1] - v[i, J1 - 1]
        v[i, J2 - 1] = 2.0 * v[i, J2] - v[i, J2 + 1]
    end

    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 < i < I2) && (J1 < j < J2)
            elseif j == 2
                urhs[i, j] -= v[i, j] * (-u[i, j+2] + 8.0 * (u[i, j+1] - u[i, j-1]) + 2.0 * u[i, j-1] - u[i, j]) / (12.0*dy) +
                abs(v[i, j]) * (u[i, j+2] - 4.0 * u[i, j+1] + 6.0 * u[i, j] - 4.0 * u[i, j-1] + 2.0 * u[i, j-1] - u[i, j]) / (4.0*dy)
                vrhs[i, j] -= v[i, j] * (-v[i, j+2] + 8.0 * (v[i, j+1] - v[i, j-1]) + 2.0 * v[i, j-1] - v[i, j]) / (12.0*dy) +
                abs(v[i, j]) * (v[i, j+2] - 4.0 * v[i, j+1] + 6.0 * v[i, j] - 4.0 * v[i, j-1] + 2.0 * v[i, j-1] - v[i, j]) / (4.0*dy)
            elseif j == MY-1
                urhs[i, j] -= v[i, j] * (-2.0 * u[i, j+1] + u[i, j] + 8.0 * (u[i, j+1] - u[i, j-1]) + u[i, j-2]) / (12.0*dy) +
                abs(v[i, j]) * (2.0 * u[i, j+1] - u[i, j] - 4.0 * u[i, j+1] + 6.0 * u[i, j] - 4.0 * u[i, j-1] + u[i, j-2]) / (4.0*dy)
                vrhs[i, j] -= v[i, j] * (-2.0 * v[i, j+1] + v[i, j] + 8.0 * (v[i, j+1] - v[i, j-1]) + v[i, j-2]) / (12.0*dy) +
                abs(v[i, j]) * (2.0 * v[i, j+1] - v[i, j] - 4.0 * v[i, j+1] + 6.0 * v[i, j] - 4.0 * v[i, j-1] + v[i, j-2]) / (4.0*dy)
            else
                urhs[i, j] -= v[i, j] * (-u[i, j+2] + 8.0 * (u[i, j+1] - u[i, j-1]) + u[i, j-2]) / (12.0*dy) +
                abs(v[i, j]) * (u[i, j+2] - 4.0 * u[i, j+1] + 6.0 * u[i, j] - 4.0 * u[i, j-1] + u[i, j-2]) / (4.0*dy)
                vrhs[i, j] -= v[i, j] * (-v[i, j+2] + 8.0 * (v[i, j+1] - v[i, j-1]) + v[i, j-2]) / (12.0*dy) +
                abs(v[i, j]) * (v[i, j+2] - 4.0 * v[i, j+1] + 6.0 * v[i, j] - 4.0 * v[i, j-1] + v[i, j-2]) / (4.0*dy)
            end
        end
    end

    for i in 2:MX-1
        for j in 2:MY-1
            if (I1 <= i <= I2) && (J1 <= j <= J2)
            else
                u[i, j] += dt * urhs[i, j]
                v[i, j] += dt * vrhs[i, j]
            end
        end
    end
end

function main()
    CFL = 0.2

    dx = 0.1
    dy = 0.1
    dt = CFL * min(dx, dy)
    nbegin = 0
    time = 0.0
    u = ones(MX, MY)
    v = zeros(MX, MY)
    p = zeros(MX, MY)

    BoundaryConditionforP(p)
    BoundaryConditionforV(u, v)

    CSV.write("logp.csv", DataFrame(p))
    CSV.write("logu.csv", DataFrame(u))
    CSV.write("logv.csv", DataFrame(v))

    println("Created files to save csv.")

    for i in 1:MAXSTEP
        nstep = i + nbegin
        time += dt
        PoissonEquation(u, v, p, dx, dy, dt)
        BoundaryConditionforP(p)
        VelocityEquation(u, v, p, dx, dy, dt)
        BoundaryConditionforV(u, v)
        println(i)
        if i % IO_FREQUENCY == 0
            CSV.write("logp.csv", DataFrame(p); append=true)
            CSV.write("logu.csv", DataFrame(u); append=true)
            CSV.write("logv.csv", DataFrame(v); append=true)
            print("Saved csv at step ", i, ".\n")
        end
    end
end

main()
