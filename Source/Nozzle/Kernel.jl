#main function to perform Chemically-coupled MOC (CHEMOC)
include("UnitProcess_INTERIOR_Cminus.jl")
include("UnitProcess_AXIS.jl");
include("Initial_line.jl");
#this module computes the kernel
function compute_kernel(IV)
        n=length(IV);
        #create empty vector to push values into it:
        V=[];
        #push the first value:
        push!(V,[IV[1]]);
        #start the i loop from 2:n
        for i = 2:n
                #create a vector of length i
                V_i=Vector{Dict{Symbol, Any}}(undef,i);
                for j= 1:i
                        if j==1
                                V_i[j]=IV[i];
                        elseif j!=i
                                # print("\nNode2 properties:\n",V_i[j-1],"\nGAS:",V_i[j-1][:gas]());
                                # print("\nNode3 properties:\n",V[i-1][j-1],"\nGAS:",V[i-1][j-1][:gas]());
                                # print("\nNode5 properties:\n",V[i-1][j],"\nGAS:",V[i-1][j][:gas]());
                                # print("\nNode5 gas:\n",V[i-1][j][:gas].X);
                                V_i[j]=compute_new_Cminus(V_i[j-1],V[i-1][j-1],V[i-1][j]);
                        else
                                V_i[j]=compute_new_axis(V_i[i-1],V[i-1][j-1]);
                        end
                end
                push!(V,V_i);
        end
        return V
end
#main:TEST()
n=10;
init=initial_line(1,n);
res_vec=compute_kernel(init);
plot(title="Kernel for $n points");
x=[];y=[];
for i =1:n
        for j=1:i
                push!(x,res_vec[i][j][:x]);
                push!(y,res_vec[i][j][:y]);
        end
end
display(scatter(x,y));
#find the mach number along the axis
M=[];x_ax=[];C=[];P=[];V=[];Temp=[];
for i=1:n
        push!(M,1/sin(res_vec[i][i][:α]));
        push!(C,res_vec[i][i][:gas].X[2]);
        push!(P,res_vec[i][i][:gas].P);
        push!(Temp,res_vec[i][i][:gas].T);
        push!(V,res_vec[i][i][:V]);
        push!(x_ax,res_vec[i][i][:x]);
end
display(plot(x_ax,M));
