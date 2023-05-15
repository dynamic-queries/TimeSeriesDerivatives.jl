abstract type AbstractDerivative end
struct TVD <: AbstractDerivative end 
abstract type FDOrder<: AbstractDerivative end 
struct IIc <: FDOrder end 
struct IVc <: FDOrder end  
struct VIc <: FDOrder end 
struct VIIIc <: FDOrder end 
struct If <: FDOrder end
struct IIf <: FDOrder end
struct IIIf <: FDOrder end 
struct IVf <: FDOrder end 
struct Vf <: FDOrder end 
struct VIf <: FDOrder end 
struct Ib <: FDOrder end 
struct IIb <: FDOrder end 
struct IIIb <: FDOrder end 

mutable struct TSDerivative
    x_raw::AbstractArray
    dt::Float64
    x::Any
    x_dot::Any

    function TSDerivative(x::AbstractArray,dt::Float64)
        new(x,dt,nothing,nothing)
    end 
end 

function get_inferables(tsd::TSDerivative)
    tsd.x, tsd.x_dot
end 

# https://en.wikipedia.org/wiki/Finite_difference_coefficient 

#----------------Central Difference--------------------#

function (tsd::TSDerivative)(::IIc)
    tsd.x = tsd.x_raw[2:end-1]
    tsd.x_dot = (0.5/tsd.dt) * (tsd.x_raw[3:end] .- tsd.x_raw[1:end-2])  
    nothing
end 

function (tsd::TSDerivative)(::IVc)
    tsd.x = tsd.x_raw[3:end-2]
    tsd.x_dot = (1/tsd.dt) *( 
                                +(1/12)*tsd.x_raw[1:end-4]
                                -(2/3)*tsd.x_raw[2:end-3]
                                +(2/3)*tsd.x_raw[4:end-1]
                                -(1/12)*tsd.x_raw[5:end]
                            )
    nothing
end 

function (tsd::TSDerivative)(::VIc)
    tsd.x = tsd.x_raw[4:end-3]
    tsd.x_dot = (1/tsd.dt) * (
                                -(1/60)*tsd.x_raw[1:end-6]
                                +(3/20)*tsd.x_raw[2:end-5]
                                -(3/4)*tsd.x_raw[3:end-4]
                                +(3/4)*tsd.x_raw[5:end-2]
                                -(3/20)*tsd.x_raw[6:end-1]
                                +(1/60)*tsd.x_raw[7:end]
                             )
    nothing
end 

function (tsd::TSDerivative)(::VIIIc)
    tsd.x = tsd.x_raw[5:end-4]
    tsd.x_dot = (1/tsd.dt) * (  
                                +(1/280)*tsd.x_raw[1:end-8]
                                -(4/105)*tsd.x_raw[2:end-7]
                                +(1/5)*tsd.x_raw[3:end-6]
                                -(4/5)*tsd.x_raw[4:end-5]
                                +(4/5)*tsd.x_raw[6:end-3]
                                -(1/5)*tsd.x_raw[7:end-2]
                                +(4/105)*tsd.x_raw[8:end-1]
                                -(1/280)*tsd.x_raw[9:end]
                             ) 

    nothing
end 

#----------------Forward Difference--------------------#

function (tsd::TSDerivative)(::If)
    tsd.x = tsd.x_raw[1:end-1]
    tsd.x_dot = (1/tsd.dt) * (
                                tsd.x_raw[2:end] - tsd.x_raw[1:end-1]
                             )
    nothing
end 

function (tsd::TSDerivative)(::IIf)
    tsd.x = tsd.x_raw[1:end-2]
    tsd.x_dot = (1/tsd.dt) * (
                                -(3/2)*tsd.x_raw[1:end-2]
                                +(2)*tsd.x_raw[2:end-1]
                                -(1/2)*tsd.x_raw[3:end]
                             )
    nothing
end

function (tsd::TSDerivative)(::IIIf)
    tsd.x = tsd.x_raw[1:end-3]
    tsd.x_dot = (1/tsd.dt) * (
                                -(11/6)*tsd.x_raw[1:end-3]
                                +(3)*tsd.x_raw[2:end-2]
                                -(3/2)*tsd.x_raw[3:end-1]
                                +(1/3)*tsd.x_raw[4:end]
                             )
    nothing
end 

function (tsd::TSDerivative)(::IVf)
    tsd.x = tsd.x_raw[1:end-4]
    tsd.x_dot = (1/tsd.dt) * (
                                -(25/12)*tsd.x_raw[1:end-4]
                                +(4)*tsd.x_raw[2:end-3]
                                -(3)*tsd.x_raw[3:end-2]
                                +(4/3)*tsd.x_raw[4:end-1]
                                -(1/4)*tsd.x_raw[5:end]
                             )
    nothing
end 

function (tsd::TSDerivative)(::Vf)
    tsd.x = tsd.x_raw[1:end-5]
    tsd.x_dot = (1/tsd.dt) * (
                                -(137/60)*tsd.x_raw[1:end-5]
                                +(5)*tsd.x_raw[2:end-4]
                                -(5)*tsd.x_raw[3:end-3]
                                +(10/3)*tsd.x_raw[4:end-2]
                                -(5/4)*tsd.x_raw[5:end-1]
                                +(1/5)*tsd.x_raw[6:end]
                             )
    nothing
end 

function (tsd::TSDerivative)(::VIf)
    tsd.x = tsd.x_raw[1:end-6]
    tsd.x_dot = (1/tsd.dt) * (
                                -(49/20)*tsd.x_raw[1:end-6]
                                +(6)*tsd.x_raw[2:end-5]
                                -(15/2)*tsd.x_raw[3:end-4]
                                +(20/3)*tsd.x_raw[4:end-3]
                                -(15/4)*tsd.x_raw[5:end-2]
                                +(6/5)*tsd.x_raw[6:end-1]
                                -(1/6)*tsd.x_raw[7:end]
                             )
    nothing    
end 

#----------------Backward Difference--------------------#

function (tsd::TSDerivative)(::Ib)
    tsd.x = tsd.x_raw[2:end]
    tsd.x_dot = (1/tsd.dt) * (
                                -1*tsd.x_raw[1:end-1]
                                +1*tsd.x_raw[2:end]
                             )
    nothing    
end 

function (tsd::TSDerivative)(::IIb)
    tsd.x = tsd.x_raw[3:end]
    tsd.x_dot = (1/tsd.dt) * (
                                +(1/2)*tsd.x_raw[1:end-2]
                                -(2)*tsd.x_raw[2:end-1]
                                +(3/2)*tsd.x_raw[3:end]
                             )
    nothing
end 

function (tsd::TSDerivative)(::IIIb)
    tsd.x = tsd.x_raw[4:end]
    tsd.x_dot = (1/tsd.dt) * (
                                -(1/3)*tsd.x_raw[1:end-3]
                                +(3/2)*tsd.x_raw[2:end-2]
                                -(3)*tsd.x_raw[3:end-1]
                                +(11/6)*tsd.x_raw[4:end]
                             )
    nothing
end 

#----------------Total variational derivative--------------------#

function (tsd::TSDerivative)(::TVD)

end