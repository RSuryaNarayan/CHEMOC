#unit process file for an axis point
#chemistry functions
include("Chemistry.jl");
include("Coefficients.jl");
#wall boundary
function wall(x)
    #returns in mm
    return (22.1852+0.71568*x-0.0010787*x^2)
end
#derivative of the wall boundary
function der_wall(x)
    #x should be in mm
    return 0.71568-0.0021574*(x)
end
#compute wall point for one iteration
function compute_new_point_drwall(node2,node3,node4,node_plus)
    #initialize α everywhere
    mach(node2);mach(node3);mach(node4);mach(node_plus);
    #compute location
    λplus=tan(node_plus[:θ]+node_plus[:α]);
    #Nlsolve
    function f!(F,x)
        F[1]=wall(x[1]*1000)-λplus*x[1]*1000-node2[:y]*1000+λplus*node2[:x]*1000;
    end
    res=nlsolve(f!,[node2[:x]],show_trace=false);
    node4[:x]=res.zero[1];
    node4[:y]=wall(node4[:x]*1000)/1000;
    node4[:θ]=atan(der_wall(node4[:x]));
    #Thermodynamics and flow properties at (4)
    Q_plus=Q(node_plus);
    T_plus=T(node4,node2,node_plus,Q_plus,1);
    P4=(T_plus-node4[:θ])/Q_plus;
    node4[:V]=node3[:V]+(node3[:gas].P-P4)/node3[:gas].density/node3[:V];
    ρ4=node3[:gas].density+(J(node3))*(node3[:x]-node4[:x])+(P4-node3[:gas].P)/(af(node3[:gas]))^2;
    C4=compute_mass_fractions(node3,node4);
    node4[:gas].DPX=ρ4,P4,C4;
    mach(node4);
end
#final march function
function march_drwall(node2,node3)
    node_plus=initialize(node2);
    node4=initialize(node2);
    #set up the convergence monitor
    n=20;iterations=collect(1:1:n);pressure=zeros(n);
    density=zeros(n);temp=zeros(n);x=zeros(n);y=zeros(n);θ=zeros(n);V=zeros(n);α=zeros(n);
    for i=1:n
        compute_new_point_drwall(node2,node3,node4,node_plus)
        pressure[i]=node4[:gas].P;
        density[i]=node4[:gas].density;
        temp[i]=node4[:gas].T;
        x[i]=node4[:x];
        y[i]=node4[:y];
        θ[i]=node4[:θ];
        V[i]=node4[:V];
        α[i]=node4[:α];
        node_plus=update(node_plus,node4)
    end
    #Convergence monitor
    # p1=plot(iterations,pressure,title="Thermodynamics vs. Iterations", xlabel="Iterations", ylabel="Pressure",linecolor="green",lw=3);
    # p2=plot(iterations,density,linecolor="red",lw=3, label="Density",xlabel="Iterations", ylabel="Density");
    # p3=plot(iterations,temp,linecolor="yellow",lw=3, label="Temperature",xlabel="Iterations", ylabel="Temperature");
    # p4=plot(iterations,x,linecolor="blue",lw=3, label="X-coordinate",xlabel="Iterations", ylabel="X");
    # p5=plot(iterations,y,linecolor="orange",lw=3, label="Y-coordinate",xlabel="Iterations", ylabel="Y");
    # p6=plot(iterations,θ,linecolor="blue",lw=3,label="θ",xlabel="Iterations",ylabel="θ");
    # p7=plot(iterations,V,linecolor="black",lw=3,label="V",xlabel="Iterations",ylabel="V");
    # p8=plot(iterations,α,linecolour="magenta",lw=3,label="α",xlabel="Iterations",ylabel="α");
    # plt1=plot(p1,p2,layout=(2,1), legend=false); display(plt1);
    # plt2=plot(p3,p4,layout=(2,1),legend=false); display(plt2);
    # plt3=plot(p5,p6,layout=(2,1),legend=false);display(plt3);
    # plt4=plot(p7,p8,layout=(2,1),legend=false);display(plt4);
    return node4
end
#test main()
#throat radius and arc radius
# r_t=0.5;r=0.5;y_c=r_t+r;
# gas2=ct.Solution("air.xml");
# gas2.TPX=1850,50000,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
# node2=Dict(:x=>60.480/1000,:y=>59.625/1000,:V=>2274.2,:θ=>30.122*pi/180,:gas=>gas2);
# gas3=ct.Solution("air.xml");
# gas3.TPX=1500,50000,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
# node3=Dict(:x=>55.943/1000,:y=>58.845/1000,:V=>2252.9,:θ=>30.752*pi/180,:gas=>gas3);
# gas=ct.Solution("air.xml");
# #dummy nodes
# gas.TP=300,101325;
# node4=Dict(:x=>1,:y=>1,:V=>1000,:θ=>10*pi/180,:gas=>gas);
# node4=march_drwall(node2,node3);
