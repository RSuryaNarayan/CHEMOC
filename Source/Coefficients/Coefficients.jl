#A subroutine that includes all the coefficients used in the unit processes
#mach angle
function mach(node)
    node[:α]=asin(af(node[:gas])/node[:V]);
end
#some predictor functions to make it easy for computation
function Q(node)
    M=node[:V]/af(node[:gas]);
    ρ=node[:gas].density;
    V=node[:V];
    return sqrt(M^2-1)/ρ/V^2;
end
#another predictor function
function S(node,c)
    θ=node[:θ];α=node[:α];y=node[:y];M=node[:V]/af(node[:gas]);
    ϕ=ψ(node[:gas]);ρ=node[:gas].density;V=node[:V];a=af(node[:gas]);
    if c==1
        return (sin(θ)/y/M-ϕ/ρ/V^2/a)*(1/cos(θ+α))
    else
        return (sin(θ)/y/M-ϕ/ρ/V^2/a)*(1/cos(θ-α))
    end
end
#yet another predictor function
function T(node4,node2,nodeplus,q,c)
    s=S(nodeplus,c);x4=node4[:x];x2=node2[:x];
    p=node2[:gas].P;θ=node2[:θ];
    if c==1
        return -s*(x4-x2)+q*p+θ
    else
        return -s*(x4-x2)+q*p-θ
    end
end
#yet another predictor function
function J(node)
    return ψ(node[:gas])/node[:V]/cos(node[:θ])/(af(node[:gas]))^2
end
#initialize the predictor
function initialize(node)
    gas_plus=ct.Solution("air.xml");
    gas_plus.DPX=node[:gas].density,node[:gas].P,node[:gas].X;
    node_plus=Dict(:x=>node[:x],:y=>node[:y],:V=>node[:V],:θ=>node[:θ],:gas=>gas_plus);
    return node_plus
end
#to apply the corrector
function update(node_plus,node1,node4)
    P4=0.5*(node1[:gas].P+node4[:gas].P);
    ρ4=0.5*(node1[:gas].density+node4[:gas].density);
    C4=0.5*(node1[:gas].X+node4[:gas].X);
    node_plus[:gas].DPX=ρ4,P4,C4;
    node_plus[:θ]=0.5*(node1[:θ]+node4[:θ]);
    node_plus[:V]=0.5*(node1[:V]+node4[:V]);
    node_plus[:y]=0.5*(node1[:y]+node4[:y]);
    mach(node_plus);
    return node_plus
end
