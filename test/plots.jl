using LodaAnomaly, LinearAlgebra, Statistics, PGFPlotsX, ADatasets
using LodaAnomaly: Loda, optimize_loda_params, criterion
using ADatasets: auc, roccurve

normaldata(n, ϵ) = cholesky([2.0 ϵ; ϵ 1.0]).L * randn(2, n)
clustered_anomalies(n,σ = 0.1) = σ .* randn(2,n) .+ (3.5,3.5)
function scattered_anomalies(n)
	r = 0.9 .+ 0.1 .* (rand(n));
	t = 2*pi*rand(n);
	X = 3 * r .* cos.(t);
	Y = 3 * r .* sin.(t);
	x = 1.5 .* (10*X-3*Y) ./ sqrt(58);
	y = (3*X+10*Y) ./sqrt(58);
	vcat(x',y')
end

n = 10000
na = Int(0.05 * n)
x = normaldata(1000, 0.9)

(ws, histograms, fs) = optimize_loda_params(x, 0.001);
c = [criterion(fs[1:i]) for i in 3:length(fs)]
plot(3:length(c)+2, log.(c), label = "normal");

a = scattered_anomalies(na)
(ws, histograms, fs) = optimize_loda_params([x a], 0.001);
cs = [criterion(fs[1:i]) for i in 3:length(fs)]
plot!(3:length(cs)+2, log.(cs), label = "scattered");

a = clustered_anomalies(na)
(ws, histograms, fs) = optimize_loda_params([x a], 0.001);
cc = [criterion(fs[1:i]) for i in 3:length(fs)]
plot!(3:length(cc)+2, log.(cc), label = "clustered")

###############################################################################
# how number of projections depend on the fraction of anomalies
###############################################################################
n = 10000
αs = vcat(0, 0.1 .* 2 .^ (1:10) ./ 1000)
truelabels = vcat(zeros(Int,1000),ones(Int,1000))
ks = map([scattered_anomalies, clustered_anomalies]) do genanom
	map(αs) do α
		na = Int(α * n)
		mean((
			x = normaldata(1000, 0.9);
			a = genanom(na);
			m = Loda([x a], 0.01);
			e = auc(roccurve(m([normaldata(1000, 0.9) scattered_anomalies(1000)]), truelabels)...);
			[length(m.ws), e]) for _ in 1:100)
	end
end
aucs = map(ts -> [t[2] for t in ts])
ks = map(ts -> [t[1] for t in ts])

p = @pgf TikzPicture(
        Axis({ymin = 0, xmode = {log}, xlabel = "fraction of anomalies", ylabel = "number of projections",
        	"legend style"= {at="(0.98,0.02)", anchor = "south east"}
        	},
            PlotInc(Table(; x = αs, y = ks[1])),
            PlotInc(Table(; x = αs, y = ks[2])),
        Legend(["scattered","clustered"])))
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/ks_alpha.pdf", p)
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/ks_alpha.tex", p)

p = @pgf TikzPicture(
        Axis({ymin = 0, xmode = {log}, xlabel = "fraction of anomalies", ylabel = "area under ROC curve",
        	"legend style"= {at="(0.98,0.02)", anchor = "south east"}
        	},
            PlotInc(Table(; x = αs, y = aucs[1])),
            PlotInc(Table(; x = αs, y = aucs[2])),
        Legend(["scattered","clustered"])))
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/aucs_alpha.pdf", p)
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/aucs_alpha.tex", p)

###############################################################################
# how number of projections depend on tau
###############################################################################
n = 10000
τs = 0.1 .* 2 .^ (5:10) ./ 1000
ks = map([scattered_anomalies, clustered_anomalies]) do genanom
	map(τs) do τ
		na = Int(0.05 * n)
		mean((
		x = normaldata(1000, 0.9);
		a = genanom(na);			
		m = Loda([x a], 0.01);
		e = auc(roccurve(m([normaldata(1000, 0.9) scattered_anomalies(1000)]), truelabels)...);
		[length(m.ws), e]) for _ in 1:100)
	end
end

aucs = map(ts -> [t[2] for t in ts])
ks = map(ts -> [t[1] for t in ts])

p = @pgf TikzPicture(
        Axis({ymin = 0, xmode = {log}, xlabel = "\$\\tau\$", ylabel = "number of projections",
        	"legend style"= {at="(0.98,0.98)", anchor = "north east"}
        	},
            PlotInc(Table(; x = τs, y = ks[1])),
            PlotInc(Table(; x = τs, y = ks[2])),
        Legend(["scattered","clustered"])))
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/ks_tau.pdf", p)
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/ks_tau.tex", p)

p = @pgf TikzPicture(
        Axis({ymin = 0, xmode = {log}, xlabel = "\$\\tau\$", ylabel = "number of projections",
        	"legend style"= {at="(0.98,0.98)", anchor = "north east"}
        	},
            PlotInc(Table(; x = τs, y = aucs[1])),
            PlotInc(Table(; x = τs, y = aucs[2])),
        Legend(["scattered","clustered"])))
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/aucs_tau.pdf", p)
PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/aucs_tau.tex", p)

###############################################################################
# Scatter plot  of problems
###############################################################################
x = normaldata(1000, 0.9)
for (f, name) in zip([scattered_anomalies, clustered_anomalies], ["scattered", "clustered"])
	a = f(100)
	p = @pgf TikzPicture(
	        Axis({ymin = -4.5, ymax = 4.5, xmin = -6, xmax = 6,
	        	"legend style"= {at="(0.98,0.02)", anchor = "south east", title = "$(name) anomalies"}
	        	},
	            PlotInc({"only marks"}, Table(; x = x[1,:], y = x[2,:])),
	            PlotInc({"only marks", mark="*"}, Table(; x = a[1,:], y = a[2,:])),
	        Legend(["normal","anomalous"])))
	PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/$(name).pdf", p)
	PGFPlotsX.pgfsave("/Users/tpevny/Documents/Habilitation/presentation/loda/$(name).tex", p)
end
