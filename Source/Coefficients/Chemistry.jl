#Chemistry terms occuring in the unit process and their computation
#PREAMBLE
#run the following comman if the mentioned packages have not been installed
# using Pkg;Pkg.add("PyCall");Pkg.add("NLsolve");Pkg.add("Plots");
#if you have not configured PYTHON into julia run the following command
# include("JuCant.jl")
#restart after configuring
using PyCall, Plots, NLsolve
@pyimport cantera as ct
#term 1:Frozen speed of sound
function af(gas)
    #raw frozen
    gamma_f =gas.cp/gas.cv;
    # gamma_f=1.2;
    a_frozen=sqrt(gamma_f*gas.P/gas.density);
    return a_frozen;
end
#species source term
function σ(gas)
    σ=gas.net_production_rates.*gas.molecular_weights
    return σ
end
#the non-homogenous ψ
function ψ(gas)
    γ_f=gas.cp/gas.cv;
    R_i=ct.gas_constant./gas.molecular_weights;
    h_i=gas.partial_molar_enthalpies./gas.molecular_weights;
    ψ_term=sum((γ_f*R_i*gas.T-(γ_f-1)*h_i).*σ(gas));
    return ψ_term
end
#species continuity equation solution
function compute_mass_fractions(node3,node4)
    ρ4=node4[:gas].density;P4=node4[:gas].P;
    W0=ρ4*node4[:V];
    del_x0=node4[:x]-node3[:x];k=del_x0/W0;
    #set the initial guess at 4 to be the same as that of 3
    C3=node3[:gas].X;
    #begin NLsolve
    function f!(F,C4)
            node4[:gas].DPX=ρ4,P4,C4;
            σ4=σ(node4[:gas]);
            for i=1:8
                F[i]=C4[i]-C3[i]-k*σ4[i];
            end
    end
    res_vec=nlsolve(f!,C3,show_trace=false);
    return res_vec.zero;
end
#List of reactions provided for reference
# for (i,r) in enumerate(gas_test.reactions())
#     print("\n")
#     print(gas_test.reaction(i-1).equation)
# end
#testing if functions work alright
# gas_test=ct.Solution("air.xml");
# gas_test.TPX=1734.6,56446.3,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
# print(ψ(gas_test),"\n\n");
# print(σ(gas_test));
