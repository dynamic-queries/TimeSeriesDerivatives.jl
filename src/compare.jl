x = nothing
dt = nothing

# Setup evaluation
tsd = TSDerivative(x,dt)
cd = [IIc, IVc, VIc, VIIIc]
fd = [If, IIf, IIIf, IVf, Vf, VIf]
bd = [Ib, IIb, IIIb]
methods = vcat(cd,fd,bd)
results = []
function eval_deriv(tsd,mode)
    tsd(mode())
    return get_inferables(tsd)[2]
end 

# Evaluate
for diffs in methods
    push!(results,eval_deriv(tsd,diffs))
end 

# Compare
## Find the shortest array.