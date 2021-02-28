#unit process for point (3) on the axis
include("Chemistry.jl");
include("Coefficients.jl");
#compute_new_point_axis
function compute_new_point_axis(node1,node3,node4,node_minus,node_plus)
    #initialize α everywhere
    mach(node1);mach(node3);mach(node4);mach(node_minus);mach(node_plus);
    #compute location
    λminus=tan(node_minus[:θ]-node_minus[:α]);
    #set y4 and θ4
    node4[:y]=0;node4[:θ]=0;
    #compute x4;
    node4[:x]=node1[:x]-node1[:y]/λminus;
    #thermodynamics and chemistry
    Q_minus=Q(node_minus);T_minus=T(node4,node1,node_minus,Q_minus,-1);
    P4=T_minus/Q_minus;
    R0=node_plus[:gas].density*node_plus[:V];A0=(af(node_plus[:gas]))^2;
    T01=R0*node3[:V]+node3[:gas].P;T02=node3[:gas].P-A0*node3[:gas].density;
    node4[:V]=(T01-P4)/R0;
    ρ4=node3[:gas].density+(J(node_plus))*(node3[:x]-node4[:x])+(P4-node3[:gas].P)/A0;
    node4[:gas].DP=ρ4,P4;
    C4=compute_mass_fractions(node3,node4);
    node4[:gas].DPX=ρ4,P4,C4;
    mach(node4);
end
#master function march_axis
#node3 is on the axis
function compute_new_axis(node1,node3)
    node4=initialize(node1);
    node_minus=initialize(node1);
    node_plus=initialize(node3);
    #set up the convergence monitor
    n=10;iterations=collect(1:1:n);pressure=zeros(n);
    density=zeros(n);temp=zeros(n);x=zeros(n);y=zeros(n);θ=zeros(n);V=zeros(n);α=zeros(n);
    for i=1:n
        compute_new_point_axis(node1,node3,node4,node_minus,node_plus)
        pressure[i]=node4[:gas].P;
        density[i]=node4[:gas].density;
        temp[i]=node4[:gas].T;
        x[i]=node4[:x];
        y[i]=node4[:y];
        θ[i]=node4[:θ];
        V[i]=node4[:V];
        α[i]=node4[:α];
        node_minus=update(node_minus,node1,node4);
        node_plus=update(node_plus,node3,node4);
    end
    if (abs(x[end]-x[end-1])<1e-5 &&
        abs(y[end]-y[end-1])<1e-5 &&
        abs(θ[end]-θ[end-1])<1e-5 &&
        abs(V[end]-V[end-1])<1e-5 &&
        abs(α[end]-α[end-1])<1e-5 &&
        abs(pressure[end]-pressure[end-1])<1e-5 &&
        abs(density[end]-density[end-1])<1e-5 &&
        abs(temp[end]-temp[end-1])<1e-5)
        print("\nConvergence Successful-Axis point\n")
    end
    #Convergence monitor
    p1=plot(iterations,pressure,title="Axis point", xlabel="Iterations", ylabel="Pressure",linecolor="green",lw=3);
    p2=plot(iterations,density,linecolor="red",lw=3, label="Density",xlabel="Iterations", ylabel="Density");
    p3=plot(iterations,temp,linecolor="yellow",lw=3, label="Temperature",xlabel="Iterations", ylabel="Temperature");
    p4=plot(iterations,x,linecolor="blue",lw=3, label="X-coordinate",xlabel="Iterations", ylabel="X");
    p5=plot(iterations,y,linecolor="orange",lw=3, label="Y-coordinate",xlabel="Iterations", ylabel="Y");
    p6=plot(iterations,θ,linecolor="blue",lw=3,label="θ",xlabel="Iterations",ylabel="θ");
    p7=plot(iterations,V,linecolor="black",lw=3,label="V",xlabel="Iterations",ylabel="V");
    p8=plot(iterations,α,linecolour="magenta",lw=3,label="α",xlabel="Iterations",ylabel="α");
    plt1=plot(p1,p2,layout=(2,1), legend=false); display(plt1);
    plt2=plot(p3,p4,layout=(2,1),legend=false); display(plt2);
    plt3=plot(p5,p6,layout=(2,1),legend=false);display(plt3);
    plt4=plot(p7,p8,layout=(2,1),legend=false);display(plt4);
    return node4;
end
#main:test()
# gas1=ct.Solution("air.xml");
# gas1.DP=0.32947,1.7025e5;
# node1=Dict(:x=>79.625/1000,:y=>1.290/1000,:V=>2306.3,:θ=>0.8866*pi/180,:gas=>gas1);
# gas3=ct.Solution("air.xml");
# gas3.DP=0.36225,1.9077e5;
# node3=Dict(:x=>76.195/1000,:y=>0/1000,:V=>2280.4,:θ=>0*pi/180,:gas=>gas3);
# gas=ct.Solution("air.xml");
# res=compute_new_axis(node1,node3);
