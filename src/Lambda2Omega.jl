module Lambda2Omega

using Interpolations

export wavelength_range
export linear_omega_range,lambda2omega

"Convenience function for defining wavelength calibration of a spectrometer."
wavelength_range(lambda0,slope,N) = StepRangeLen(lambda0,slope,N)

const cc=3e8*1e9/1e12 # speed of light in nm/ps

raw_omega_range(wl::AbstractArray) = 2*pi*cc./wl

"The linearly spaced, increasing omega range corresponding to a wavelength range."
function linear_omega_range(wl::AbstractArray)
    incr = (wl[2]-wl[1]>0.0)
    omega0 = incr ? 2*pi*cc/wl[end] : 2*pi*cc/wl[1]
    omega1 = incr ? 2*pi*cc/wl[1] : 2*pi*cc/wl[end]
    range(omega0,stop=omega1,length=length(wl))
end

const gl = Gridded(Linear())
const gc = Gridded(Constant())

"""Interpolates data as a function of wavelength to an evenly spaced, increasing frequency
axis. For multidimensional arrays, the first dimensional is assumed to be
wavelength and the others spatial."""
function lambda2omega(wl::AbstractArray,A::Array{Float64,1})
    if wl[2]-wl[1]>0
        omega=raw_omega_range(wl)[end:-1:1]
        i=interpolate((omega,),A[end:-1:1],gl)
    else
        omega=raw_omega_range(wl)
        i=interpolate((omega,),A,gl)
    end
    out=i(linear_omega_range(wl))
end

function lambda2omega(omega::AbstractArray,omegalin::AbstractArray,A::Array{Float64,1})
    i=interpolate((omega,),A,gl)
    out=i(omegalin)
end

function lambda2omega(wl::AbstractArray,A::Array{Float64,2})
    N=size(A)[2]
    if wl[2]-wl[1]>0
        omega=raw_omega_range(wl)[end:-1:1]
        i=interpolate((omega,1:N),A[end:-1:1,:],(gl,Gridded(Constant())))
    else
        omega=raw_omega_range(wl)
        i=interpolate((omega,1:N),A,(gl,Gridded(Constant())))
    end
    out=i[linear_omega_range(wl),1:N]
end

function lambda2omega(wl::AbstractArray,A::Array{Float64,3})
    N=size(A)[2]
    M=size(A)[3]
    if wl[2]-wl[1]>0
        omega=raw_omega_range(wl)[end:-1:1]
        i=interpolate((omega,1:N,1:M),A[end:-1:1,:,:],(gl,Gridded(Constant()),Gridded(Constant())))
    else
        omega=raw_omega_range(wl)
        i=interpolate((omega,1:N,1:M),A,(gl,Gridded(Constant()),Gridded(Constant())))
    end
    out=i[linear_omega_range(wl),1:N,1:M]
end

function lambda2omega(wl::AbstractArray, A::Array{UInt16})
    lambda2omega(wl,Array{Float64}(A))
end

end
