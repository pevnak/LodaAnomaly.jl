struct EqualHistogram{T}
	p :: Vector{T}
	mi::T
	mx::T
	δ::T
end

function EqualHistogram(x)
	mi, mx = minimum(x), maximum(x)
	mx - mi == 0 && @error("The span of histogram should be greater than zero")

	x = (x .- mi) ./ (mx - mi)

	k = histoptimal(x)[2]
	δ = 1 / k
	p = bindata(x, δ) ./ (length(x)*(mx - mi))
	EqualHistogram(p, mi, mx, δ)
end

function (m::EqualHistogram)(x)
	ii = floor.(Int, (x .- m.mi) ./ ((m.mx - m.mi) * m.δ))
	o = similar(x)
	fill!(o, log(10*eps(eltype(x))))
	mask = (ii .>0 ) .& (ii .< length(m.p))
	o[mask] .= -log.( m.p[ii[mask]] .+ 10*eps(eltype(x)))
	o
end
