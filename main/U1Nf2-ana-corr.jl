using Revise
import Pkg
Pkg.activate(".")
using LFTs

Pkg.develop(path="/home/david/.julia/personal-utils/Utils/")
using Utils

using ADerrors

filepath1 = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/trash/U1MassCheck/2023-01-20T15:25:53.056/measurements/corr_pion.txt"

filepath2 = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/trash/U1MassCheck/2023-01-20T15:25:53.056/measurements/corr_pcac.txt"

filepath1 = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/trash/U1CriticalMassCheck/2023-01-25T16:59:50.490/measurements/corr_pion.txt"
filepath2 = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTs.jl/results/trash/U1CriticalMassCheck/2023-01-25T16:59:50.490/measurements/corr_pcac.txt"

using DelimitedFiles

# abstract type AbstractFunction <: UtilsFunc end
# abstract type AbstractCorrelationFunction end

# mutable struct U1PionCorrelationFunction <: AbstractCorrelationFunction
#     pathfile
#     xdata
#     ydata::Vector{uwreal}
#     T
#     partner
#     function U1PionCorrelationFunction(; pathfile, partner::DataType = Any, burnout::Int64 = 1)
#         corr = read(partner, pathfile)[burnout:end, :]
#         T = length(corr[1,:])
        

#         return new(wdir,) 
#     end
# end
# nparameters(s::U1PionCorrelationFunction) = 2
# (s::U1PionCorrelationFunction)(x, p) = p[1] * (exp(-p[2]*x) + exp(-p[2]*(s.T-x)))


function extract_symmetrize(corr_file::String, ID::String; burnout::Int64=1) #{{{
	# Load correlation data
    corr = readdlm(corr_file, ',', Float64)[burnout:end,:]
	T = length(corr[1,:]);
	nconf = length(corr[:,1]);

	println("Size: ", size(corr))
	println("Temporal length: ", T)
	
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);

	for i in eachindex(ydata)
		ydata[i] = uwreal( (corr[:,i] .+ corr[:,1+(T-i+1)%T])/2, ID)
		uwerr(ydata[i])
	end

	return xdata, ydata, T
end #}}}


function plot_correlator(corr_file::String, ID::String; burnout::Int64=1) #{{{
	# Load correlation data
    corr = readdlm(corr_file, ',', Float64)[burnout:end,:]
	T = length(corr[1,:]);

	println("Size: ", size(corr))
	println("Temporal length: ", T)
	
	xdata = range(0, length=T)
	ydata = Vector{uwreal}(undef, T);

	for i in eachindex(ydata)
		ydata[i] = uwreal(corr[:,i], ID)
		uwerr(ydata[i])
	end

	pl = plot(xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, yscale=:identity)
	display(pl)

    return xdata, ydata, T

end #}}}

function tmin_fit(f::UtilsFunc , xdata, ydata::Vector{uwreal}, tmin::Int64, prms0::Vector{Float64}; Tmax::Int64=length(ydata)-1, plot_tmin_fit::Bool=false, yscale::Symbol=:identity) #{{{

	if (Tmax > xdata[end][1])
		throw("Tmax value greater than greatest value in xdata")
	end

	pos_tmin = findfirst(x->x[1]==tmin, xdata)
	pos_Tmax = findfirst(x->x[1]==Tmax, xdata)
	fit_region = pos_tmin:pos_Tmax

	# Fit
	fitp, cse, cs = fit_data(f, xdata[fit_region], ydata[fit_region], prms0)
	uwerr.(fitp) # compute errors

	# Plot fit
	if(plot_tmin_fit==true)
		if(size(xdata[1],1)==1)
			pl = plot_fit(xdata, ydata, f, fitp)
			plot!(pl, xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, title="tmin = "*string(tmin), yscale=yscale)
			# plot!(ylim=(9000,10000))
			display(pl)
		else
			res = f.(xdata, [fitp for i in 1:length(xdata)])
			uwerr.(res)
			pl_x = (hcat(xdata...) |> permutedims)[:,1]
			pl = plot( pl_x, value.(res), ribbons=err.(res), reuse=false, title="tmin="*string(tmin)*", 1")
			plot!(pl, pl_x, value.(ydata), yerr=err.(ydata), seriestype=:scatter)
			display(pl)
		end
	end
	
	return fitp, cse, cs

end #}}}

function tmin_loop(f::UtilsFunc , xdata, ydata, tmin::Int64, tmax::Int64, prms0::Vector{Float64}; update_prms::Bool=true, Tmax::Int64=length(ydata)-1, plot_tmin_fit::Bool=false, plot_column::Int64=0) #{{{

	fitps = [] # tmin and parameters to be returned
	tminvalues = tmin:tmax

	prms = prms0

	for itmin in tminvalues

		# Fit
		fitp, cse, cs = tmin_fit(f, xdata, ydata, itmin, prms, Tmax=Tmax, plot_tmin_fit=plot_tmin_fit, yscale=:log10)
		display(cs)

		if(itmin == tmin)
			fitps = hcat(itmin, permutedims( fitp ) )
		else
			aux_fitps = hcat(itmin, permutedims( fitp ))
			fitps = vcat(fitps, aux_fitps)
		end

		if(update_prms)
			prms = value.(fitp)
		end
	end

	if(plot_column > 0)
		pl = plot(fitps[:,1], value.(fitps[:,plot_column]), yerr=err.(fitps[:,plot_column]), seriestype=:scatter)
		display(pl)
	end

	return fitps
end #}}}


burnout = 1000
strID = "test"

corrs = read(U1PionCorrelator, filepath)

xdata, ydata, T = extract_symmetrize(filepath, strID, burnout=burnout)

f = SymCorrelator(T, 1, 0)

fitp, cse, cs = tmin_fit(f, xdata, ydata, 3, [0.1, 0.2])

tmin_loop

fitps = tmin_loop(f, xdata, ydata, 1,12, [0.1, 0.2], plot_column = 3)




ws = CorrelatorAnalysis(filepath1, U1PionCorrelator, burnout = 500, ensemble_ID = "test3")
uwrealsym(ws)
ws.tmax = 12

tmin_loop(ws, [1.0, 1.0])

plot(ws.histories, seriestype=:scatter)

plot(ws, seriestype=:scatter)


ws2 = CorrelatorAnalysis(filepath2, U1PCACCorrelator, burnout = 500, ensemble_ID = "test3")
uwrealsym(ws2)
ws2.tmax = 12

derivate_sym_correlator!(ws2)

tmin_loop(ws, [1.0, 1.0])

mPCAC = 1/2 * ws2.ydata ./ ws.ydata[2:end-1]

uwerr.(mPCAC)

plot(1:length(mPCAC), mPCAC, seriestype=:scatter)
