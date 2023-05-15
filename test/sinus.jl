using Test
using TimeSeriesDerivatives
using Plots
using Measures

# Artificial data
x = 0.0:1e-3:2*pi
y = sin.(x)
dx = x[2]-x[1]

# Compute derivatives

## Central differences
tsd = TSDerivative(y,dx)

tsd(IIc())
c2x,c2v = get_inferables(tsd)
tsd(IVc())
c4x,c4v = get_inferables(tsd)
tsd(VIc())
c6x,c6v = get_inferables(tsd)
tsd(VIIIc())
c8x,c8v = get_inferables(tsd)

## Forward differences
tsd(If())
f1x,f1v = get_inferables(tsd)
tsd(IIf())
f2x,f2v = get_inferables(tsd)
tsd(IIIf())
f3x,f3v = get_inferables(tsd)
tsd(IVf())
f4x,f4v = get_inferables(tsd)
tsd(Vf())
f5x,f5v = get_inferables(tsd)
tsd(VIf())
f6x,f6v = get_inferables(tsd)

## Backward differences
tsd(Ib())
b1x,b1v = get_inferables(tsd)
tsd(IIb())
b2x,b2v = get_inferables(tsd)
tsd(IIIb())
b3x,b3v = get_inferables(tsd)

# Visualize the derivatives
f2 = plot(x,cos.(x),title="Analytical Derivative")
f3 = plot(c2v, title="II")
f4 = plot(c4v, title="IV")
f5 = plot(c6v, title="VI")
f6 = plot(c8v, title="VIII")
plot(f2,f3,f4,f5,f6,layout=(3,2),size=(750,750),margin=10mm)
savefig("central.png")

# Order of convergence 
analytic = cos.(x)
e1 = c2v .- analytic[2:end-1]
e2 = c4v .- analytic[3:end-2]
e3 = c6v .- analytic[4:end-3]
e4 = c8v .- analytic[5:end-4]

g1 = plot(e1,title="II",size=(250,250))
g2 = plot(e2,title="IV",size=(250,250))
g3 = plot(e3,title="VI",size=(250,250))
g4 = plot(e4,title="VIII",size=(250,250))
plot(g1,g2,g3,g4,layout=(2,2),size=(1000,1000),margin=30mm)
savefig("error_central.png")