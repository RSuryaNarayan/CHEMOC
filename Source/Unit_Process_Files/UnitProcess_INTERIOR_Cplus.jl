#Unit Process for an interior point:Along C+ characteristic
#angles in radians for immediate parse into trig functions
#for chemistry functions
include("Initial_line.jl")
include("Chemistry.jl")
include("Coefficients.jl")
#subroutine to find properties @ 1
function interpolate_at_1(node2,node3,node5,node4,node1)
    #check if interpolation is possible
    Δy=abs(node3[:y]-node5[:y]);
    Δx=abs(node3[:x]-node5[:x]);
    if (Δy>1e-5 || Δx>1e-5)
        #iterate till convergence
        for i in 1:10
            λminus=tan(node1[:θ]-node1[:α]);
            if (Δx<1e-5)
                node1[:x]=node3[:x];
                node1[:y]=node4[:y]-λminus*(node4[:x]-node1[:x]);
            else
                λ53=(node5[:y]-node5[:x])/(node5[:x]-node3[:x]);
                coeff_mat=[1 -λ53;1 -λminus];
                const_mat=[node3[:y]-λ53*node3[:x];node4[:y]-λminus*node4[:x]];
                p=coeff_mat\const_mat;
                node1[:y]=p[1];node1[:x]=p[2];
            end
            inter=sqrt((node1[:x]-node3[:x])^2+(node1[:y]-node3[:y])^2)/sqrt((node5[:y]-node3[:y])^2+(node5[:x]-node3[:x])^2);
            node1[:θ]=node3[:θ]+(node5[:θ]-node3[:θ])*inter;
            P1=node3[:gas].P+(node5[:gas].P-node3[:gas].P)*inter;
            C1=node3[:gas].X+(node5[:gas].X-node3[:gas].X)*inter;
            ρ1=node3[:gas].density+(node5[:gas].density-node3[:gas].density)*inter;
            node1[:gas].DPX=ρ1,P1,C1;
            node1[:V]=node3[:V]+(node5[:V]-node3[:V])*inter;
            mach(node1);
        end
    else
        #no interpolation required
        node1=initialize(node3);
    end
end
#base subroutine to perform one iteration of the predictor
function compute_new_point_interior(node2,node3,node5,node4,node1,nodeplus,nodeminus,node0,c)
    mach(node2);mach(node3);mach(node5);mach(nodeplus);mach(node4);mach(node1);mach(nodeminus);
    #compute point 4 location
    λ0=tan(node0[:θ]);
    λplus=tan(nodeplus[:θ]+nodeplus[:α]);
    coeff_matrix=[1 -λ0;1 -λplus];
    const_matrix=[node3[:y]-λ0*node3[:x];node2[:y]-λplus*node2[:x]];
    p=coeff_matrix\const_matrix;
    node4[:y]=p[1];node4[:x]=p[2];
    interpolate_at_1(node2,node3,node5,node4,node1);
    #calculate p4, θ4
    Q_plus=Q(nodeplus);T_plus=T(node4,node2,nodeplus,Q_plus,1);
    if c==1
        Q_minus=Q(node1);
        T_minus=T(node4,node1,node1,Q_minus,0);
    else
        Q_minus=Q(nodeminus);
        T_minus=T(node4,node1,nodeminus,Q_minus,0);
    end
    coeff_mat=[Q_plus 1;Q_minus -1];const_mat=[T_plus;T_minus];res=coeff_mat\const_mat;
    P4=res[1];node4[:θ]=res[2];
    #calculate V4:
    R0=node0[:gas].density*node0[:V];A0=(af(node0[:gas]))^2;
    T01=R0*node3[:V]+node3[:gas].P;T02=node3[:gas].P-A0*node3[:gas].density;
    node4[:V]=(T01-P4)/R0;
    #calculate ρ4:
    ρ4=node3[:gas].density+(J(node0))*(node3[:x]-node4[:x])+(P4-node3[:gas].P)/A0;
    #mass fractions
    node4[:gas].DP=ρ4,P4;
    C4=compute_mass_fractions(node3,node4);
    node4[:gas].DPX=ρ4,P4,C4;
    mach(node4);
    return node1;
end
#master subroutine to reach the point 4
function compute_new_Cplus(node2,node3,node5)
    #initialize the solution by setting node_plus=node_2 for first iteration;
    node1=initialize(node2);node_plus=initialize(node2);
    node4=initialize(node2);node_minus=initialize(node2);
    node0=initialize(node3);
    n=20;
    iterations=collect(1:1:n);pressure=zeros(n);
    density=zeros(n);temp=zeros(n);x=zeros(n);y=zeros(n);θ=zeros(n);V=zeros(n);α=zeros(n);
    #initial iteration: pass the values at 1 and 2
    for i=1:n
        node1=compute_new_point_interior(node2,node3,node5,node4,node1,node_plus,node_minus,node0,i);
        pressure[i]=node4[:gas].P;
        density[i]=node4[:gas].density;
        temp[i]=node4[:gas].T;
        x[i]=node4[:x];
        y[i]=node4[:y];
        θ[i]=node4[:θ];
        V[i]=node4[:V];
        α[i]=node4[:α];
        node_plus=update(node_plus,node2,node4);
        node_minus=update(node_minus,node1,node4);
    end
    #Convergence monitor
    p1=plot(iterations,pressure,title="Thermodynamics vs. Iterations", xlabel="Iterations", ylabel="Pressure",linecolor="green",lw=3);
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
#PROGRAM: MAIN():test():uncomment this to test the program:
gas2=ct.Solution("air.xml");
gas2.TPX=1500,56446.3,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
node2=Dict(:x=>0.1,:y=>0.2,:V=>1.1*af(gas2),:θ=>0.01*pi/180,:gas=>gas2);
gas3=ct.Solution("air.xml");
gas3.TPX=1000,56446.3,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
node3=Dict(:x=>0.1,:y=>0.3,:V=>900.00,:θ=>0.01*pi/180,:gas=>gas3);
gas5=ct.Solution("air.xml");
gas5.TPX=1000,56446.3,[0,0.49342,0,0.00402,0.00633,0.00007,0.49615,0.00001];
node5=Dict(:x=>0.1,:y=>0.4,:V=>900.00,:θ=>0.1*pi/180,:gas=>gas5);
gas=ct.Solution("air.xml");
gas.TP=300,101325;
node1=Dict(:x=>1,:y=>1,:V=>1000,:θ=>10*pi/180,:gas=>gas);
node4=Dict(:x=>1,:y=>1,:V=>1000,:θ=>10*pi/180,:gas=>gas);
# function call definiton
# initial=initial_line(1,0.1);
node4=compute_new_Cplus(node2,node3,node5);
