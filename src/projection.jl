
struct Loda{V,E}
	ws::Vector{V}
	histograms::Vector{E}
end


criterion(fs) = sum(abs.(fs[end] .- fs[end - 1])) / sum(abs.(fs[2] .- fs[1]))

function Loda(x, τ::AbstractFloat = 0.01, sparse::Bool = false)
	sparse == true && @error("Not implemented")
	Loda(optimize_loda_params(x, τ)[1:2]...)
end

function optimize_loda_params(x,τ)
	d = size(x,1)
	ws = []
	fs = []
	histograms = []
	k = 1
	cumf = zeros(size(x,2))
	while true
		w = randn(1, d)
		p = w * x
		h = EqualHistogram(p)
		cumf .+= h(p)[:]
		push!(fs, cumf ./ k)
		push!(ws, w)
		push!(histograms, h)
		
		length(fs) > 2 && criterion(fs) < τ && break
		k += 1 
	end
	(ws, histograms, fs)
end


function (m::Loda)(x)
	o = similar(x, size(x,2))
	o .= 0
	for i in 1:length(m.ws)
		h, w = m.histograms[i], m.ws[i]
		o .-= h(w * x)[:]
	end
	o ./= length(m.ws)
end