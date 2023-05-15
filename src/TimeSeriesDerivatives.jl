module TimeSeriesDerivatives
    include("fd_stencils.jl")
    export TSDerivative, get_inferables
    export TVD
    export IIc,IVc,VIc,VIIIc
    export If,IIf,IIIf,IVf,Vf,VIf
    export Ib,IIb,IIIb

    # include("compare.jl")
    export compare
end