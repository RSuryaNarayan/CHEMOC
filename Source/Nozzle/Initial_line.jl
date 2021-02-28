#initialize the solution for the method of characteristics
#initial line at the throad is vertical with compositions and pressures taken from CEA
include("Chemistry.jl")
function initial_line(throat_radius,n)
    spacing=throat_radius/(n-1);
    initial_line_vec=Vector{Dict{Symbol, Any}}(undef,n);
    gas=ct.Solution("air.xml");
    gas.TPX=2700,34.730e5,[0.0,0.410000,0.04,0.00,0.00,0.0,0.58,0.0];
    y = collect(0:2.5:17.5);y=y/1000;
    # initialize the nodes with mach number close to 1
    for i in 1:n
        initial_line_vec[i]=Dict(:x=>0,:y=>((i-1))*spacing,:V=>1.1*af(gas),:θ=>0*pi/180,:gas=>gas);
         # initial_line_vec[i]=Dict(:x=>4.365e-3-14.832*y[i]^2,:y=>y[i],:V=>1.05*af(gas),:θ=>0*pi/180,:gas=>gas);
    end
    # initial_line_vec[1]=Dict(:x=>0.0,:y=>,:V=>1.1*af(gas),:θ=>0*pi/180,:gas=>gas);
    return initial_line_vec;
end
