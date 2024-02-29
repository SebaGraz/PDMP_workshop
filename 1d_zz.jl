function poisson_time(a, b, u)
    if b > 0
        if a < 0
            return sqrt(-log(u)*2.0/b) - a/b
        else # a[i]>0
            return sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        end
    elseif b == 0
        if a > 0
            return -log(u)/a
        else # a[i] <= 0
            return Inf
        end
    else # b[i] < 0
        if a <= 0
            return Inf
        elseif -log(u) <= -a^2/b + a^2/(2*b)
            return -sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        else
            return Inf
        end
    end
end



function pdmp(x, v, μ, σ2, T)
    t = 0.0
    trace = [(t, x, v)]
    tau =  poisson_time((x-μ)*v/σ2, v^2/σ2, rand())
    while t < T
        x += tau*v
        v *= -1  
        t += tau
        push!(trace, (t, x, v))
        tau = poisson_time((x-μ)*v/σ2, v^2/σ2, rand())
    end
    return trace
end


x, v = 0.0, 1
μ = 5.0
σ2 = 2.0
T = 100.0
xx = pdmp(x, v, μ, σ2, T)
using Plots
plot(getindex.(xx,1), getindex.(xx,2)) 
