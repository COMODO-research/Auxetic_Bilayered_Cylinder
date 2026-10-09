using Auxetic_Bilayered_Cylinder
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.Statistics
using Comodo.LinearAlgebra
using AbaqusTools
using FileIO
using Printf 

# Visualization parameters
GLMakie.closeall()

saveDir = joinpath(Auxetic_Bilayered_Cylinder_dir(),"assets","temp")
if !isdir(saveDir)
    mkdir(saveDir)      
end

fileName_inp = joinpath(saveDir, "temp.inp")

markersize1  = 10
markersize2  = 10
linewidth1   = 5
linewidth2   = 2
strokewidth1 = 2
strokewidth2 = 1.0
strokewidth3 = 0.5

function ellipse_segment(r1, r2, xs, ys, pointSpacing; np_def = 1000)
    V = [Point{3,Float64}(r1*cos(t)-r1+xs, r2*sin(t)+ys, 0.0) for t in range(0.5*pi, 0.0, np_def)]
    V = evenly_space(V, pointSpacing; close_loop = false, spline_order = 4)
    return V
end

function mirrormesh(F,V::Vector{Point{3,T}}; dir=(1,0,0)) where T<:Real
    if !iseven(sum(dir))
        F = invert_faces(F)
    end
    m = ones(T,3)
    for i in eachindex(dir)
        if dir[i]==1
            m[i] = -m[i]
        end
    end
    V = [Point{3,T}(m[1]*v[1],m[2]*v[2],m[3]*v[3]) for v in V]
    return F,V
end

function copy_xyz(F,V,C,nRep,s; flip = false)
    FT = deepcopy(F)
    VT = deepcopy(V)
    CT = deepcopy(C)
    indexShift = length(V)
    n = length(V)
    if flip
        Ff = invert_faces(F)
        Vf = [Point{3,Float64}(v[2],v[1],v[3]) for v in V]
    end
    for i = 0:1:nRep[1]-1
        for j = 0:1:nRep[2]-1
            if i==0 && j==0 

            else
                vs = Point{3,Float64}(s[1]*i,s[2]*j,0.0)                
                if flip && (iseven(i+j-1)) 
                    append!(FT,[f .+ indexShift for f in Ff])
                    append!(VT,Vf .+ vs)                
                else
                    append!(FT,[f .+ indexShift for f in F])
                    append!(VT,V .+ vs)                
                end
                append!(CT,C)
                indexShift += n
            end            
        end
    end
    return FT, VT, CT
end

element_formulation = :quadratic 

## Control parameters
nRep = (1,5) # Number of repetitions in each direction
L = π./nRep[2] # Unit cell length
A = L/2.0 # Unit cell half lenght
r1 = 0.2*L # Major radius
r2 = 0.1*L # Minor radius
shellThickness = 1.0/20.0 # Shell layer thickness
pointSpacing = L/15.0 # General point spacing for the nodes
pointSpacingSide = pointSpacing .* 1.5 # Point spacing for the sides (e.g. larger to use less nodes)
axialLength = L * nRep[1] # Longitudinal lenght

nThick = ceil(Int,shellThickness/pointSpacing) + 1 # Number of nodes across shell thickness
tolLevel = pointSpacing/100.0 # Tolerance level

tetrahedral_volume_factor = 10.0

## Derived parameters
P_front = Point{3,Float64}( axialLength/2.0, 0.0, 0.0) # Mid edge
P_back  = Point{3,Float64}(-axialLength/2.0, 0.0, 0.0) # Mid edge
P_ref = Point{3,Float64}( axialLength, 0.0, 0.0) # An external reference point 

## Define ellipse curve segment
np_def = 500 # Number of points to define ellipse
dl = A-r1-r2
xs = A/2.0-r2-(A - r1 - r2)
ys = -A/2.0
d  = 0.25* 2.0 * pi * sqrt( (r1^2 + r2^2) / 2.0)
np = ceil(Int64, d/pointSpacing)
Ve = ellipse_segment(r1, r2, xs, ys, pointSpacing; np_def = np_def)

# Create repeated unit with ellipse segments on the corners
np = ceil(Int,dl/pointSpacing)+1
Vl = range(Ve[end], Ve[end] + Point{3,Float64}(dl, 0.0, 0.0), np)
Vel = [Ve;Vl[2:end]]
Vc = deepcopy(Vel)
for (i,a) in enumerate(range(0.5*π,1.5*π,3))
     Q = RotXYZ(0.0,0.0,a)    
     if i<3
        append!(Vc,[Point{3, Float64}(Q*v) for v ∈ Vel[2:end]])        
     else
        append!(Vc,[Point{3, Float64}(Q*v) for v ∈ Vel[2:end-1]])        
     end
end
v_shift = Point{3,Float64}(-A/2.0,-A/2.0,0.0)

np = ceil(Int,r2/pointSpacing)+1
Vl2 = range(Ve[1], Ve[1]- Point{3,Float64}(0.0, r2, 0.0), np)
np = ceil(Int,r1/pointSpacing)+1
Vl3 = range(Vl2[end], Vl2[end] + Point{3,Float64}(r1, 0.0, 0.0), np)

Vc2 = [reverse(Ve); Vl2[2:end]; Vl3[2:end];]

# Triangulate the lattice (not ellise segments)
Vcd1 = Vc .+ v_shift
Vcd1p = deepcopy(Vcd1)
Ft1,Vt1,_ = regiontrimesh((Vcd1,),([1],),(pointSpacing))
invert_faces!(Ft1)

# Triangulate the ellipse segment corners
Vcd2 = Vc2  .+ v_shift
Vcd2p = deepcopy(Vcd2)
Ft2,Vt2,_ = regiontrimesh((Vcd2,),([1],),(pointSpacing))
invert_faces!(Ft2)

# Copy/rotate ellipse segment regions
Ft2_copy = deepcopy(Ft2)
Vt2_copy = deepcopy(Vt2)
for (i,a) in enumerate(range(0.5*π,1.5*π,3))
    Q = RotXYZ(0.0,0.0,a)    
    if i<3
       append!(Vt2, [Point{3, Float64}(Q*(v - v_shift)) + v_shift for v ∈ Vt2_copy])        
    else
       append!(Vt2, [Point{3, Float64}(Q*(v - v_shift)) + v_shift for v ∈ Vt2_copy])        
    end
    append!(Ft2, [f .+ i*length(Vt2_copy) for f in Ft2_copy])
end

## Now copy features to create 2x2 square
Ft12, Vt12 = mirrormesh(Ft1, Vt1; dir=(1,0,0))
Ft13, Vt13 = mirrormesh(Ft1, Vt1; dir=(1,1,0))
Ft14, Vt14 = mirrormesh(Ft1, Vt1; dir=(0,1,0))


Ft22, Vt22 = mirrormesh(Ft2, Vt2; dir=(1,0,0))
Ft23, Vt23 = mirrormesh(Ft2, Vt2; dir=(1,1,0))
Ft24, Vt24 = mirrormesh(Ft2, Vt2; dir=(0,1,0))

F_unit,V_unit,C_unit = joingeom(Ft1, Vt1, Ft12, Vt12, Ft13, Vt13, Ft14, Vt14, Ft2, Vt2, Ft22, Vt22, Ft23, Vt23, Ft24, Vt24)
C_unit[C_unit.>4] .= 5
F_unit,V_unit = mergevertices(F_unit,V_unit)

## Copy 2x2 unit into strip
FT, VT, CT = copy_xyz(F_unit, V_unit, C_unit, nRep, (L,L); flip = false)
FT, VT = mergevertices(FT, VT)

###############################################################################
# Visualise strip

## Visualise repeatable unit
fig1 = Figure(size = (1200,800), fontsize = 20)
ax1 = Axis(fig1[1, 1], aspect = DataAspect(), title="Curves for 1/8 of periodic unit")
lines!(ax1, Vcd1p[[1:length(Vcd1p); 1]], color=:red, linewidth=3.0)
lines!(ax1, Vcd2p, color=:blue, linewidth=5.0)

ax2 = Axis(fig1[1, 2], aspect = DataAspect(), title="Surface mesh for periodic units")
hp1 = meshplot!(ax2, Ft1, Vt1, strokewidth=strokewidth1, color=:steelblue1, shading=false)
hp2 = meshplot!(ax2, Ft2, Vt2, strokewidth=strokewidth1, color=:lightgray, shading=false)

FTp, VTp = separate_vertices(FT, VT)
CTp = simplex2vertexdata(FTp, CT) # Convert face color data to vertex color data 

ax3 = Axis(fig1[2, 1], aspect = DataAspect(), title="Repeated periodic units")
hp3 = meshplot!(ax3, FTp, VTp, strokewidth=0.0, color=CTp, colormap=Makie.Categorical(Makie.Reverse(:Spectral)), shading=false)

# Legend(fig1[1, 3], [hp22, hp32], ["Tetrahedral core", "Pentahedral layer"])
Colorbar(fig1[1, 3], hp3)
screen = display(GLMakie.Screen(), fig1)
###############################################################################

# Morph strip to form half cylinder
t = [pi*((v[2]+L/2)./(nRep[2]*L)) for v in VT] # Angle "theta" derived from y-coordinate
VT = [Point{3,Float64}(VT[i][1], cos(t[i]), sin(t[i])) for i in eachindex(VT)]

## Thicken into pentahedral elements 
ETb = boundaryedges(FT) 

# Shift so corner is at origin
v_shift = Point{3,Float64}(L/2,0.0,0.0)
VT .+= v_shift

F_extrude,V_extrude,_ = remove_unused_vertices(FT[CT.<5],VT)
B = [v[3] < tolLevel for v in V_extrude]
indLow = findall(B)
N = vertexnormal(F_extrude,V_extrude)

for (i,n) in enumerate(N[indLow])
    if dot(n,Vec{3,Float64}(0.0,1.0,0.0))>0.0
        N[indLow[i]] = Vec{3,Float64}(0.0,1.0,0.0)
    else
        N[indLow[i]] = Vec{3,Float64}(0.0,-1.0,0.0)
    end
end

E_penta, V_penta = extrudefaces(F_extrude,V_extrude; extent=shellThickness, direction=:positive, num_steps=nThick, N = N)
if element_formulation == :quadratic
    E_penta, V_penta = penta6_penta15(E_penta, V_penta)
end

## Construct tetrahedral interior
# Find front-arc and bottom-front boundary edges 
edges_x = Vector{LineFace{Int}}()
edges_y = Vector{LineFace{Int}}()
for (i,e) in enumerate(ETb) # For all edges        
    if all([v[1] for v in VT[e]] .< tolLevel)
        push!(edges_x,e)
    end
    if all([v[2] for v in VT[e]] .< -1.0+tolLevel)
        push!(edges_y,e)
    end
end

# Construct curves
ind_curve_x = edges2curve(edges_x)
ind_curve_y = edges2curve(edges_y)

V_curve_x0_sparse = [VT[ind_curve_x[length(ind_curve_x)]], Point{3,Float64}(0.0,0.0,0.0), VT[ind_curve_x[1]]]
V_curve_x0 = evenly_space(V_curve_x0_sparse, pointSpacing; close_loop = false, spline_order = 2, must_points = [1,2,3])      

# Meshing front surface
Vf_region_curve = [VT[ind_curve_x]; V_curve_x0[2:end-1]]
Vf_region_curve = [Point{3,Float64}(v[2],v[3],v[1]) for v in Vf_region_curve] # Change to XY space for meshing
Ff,Vf,Cf = regiontrimesh((Vf_region_curve,),([1],),(pointSpacingSide))
Vf = [Point{3,Float64}(v[3],v[1],v[2]) for v in Vf] # Change to XY space for meshing

# ###############################################################################
# # Visualise strip

# ## Visualise repeatable unit
# fig1 = Figure(size = (1200,800), fontsize = 20)
# ax1 = AxisGeom(fig1[1, 1], title="Curves")
# lines!(ax1, VT[ind_curve_x], color=:red, linewidth=2.0)
# lines!(ax1, V_curve_x0[2:end-1], color=:blue, linewidth=2.0)
# scatter!(V_curve_x0[2:end-1],color=:black, markersize=15)
# hp1 = meshplot!(ax1, Ff,Vf, strokewidth=0.0, color=:white)
# screen = display(GLMakie.Screen(), fig1)

# ###############################################################################

# Create back by copying front
Fb = deepcopy(Ff)
Vb = deepcopy(Vf)
v_shift = Point{3,Float64}(axialLength,0.0,0.0)
Vb .+= v_shift
invert_faces!(Ff)

# Construct curves for bottom surface 
Vc1 = reverse(V_curve_x0)
Vc2 = VT[ind_curve_y[2:end-1]]
v_shift = Point{3,Float64}(axialLength,0.0,0.0)
Vc3 = deepcopy(V_curve_x0)
Vc3 .+= v_shift
v_shift = Point{3,Float64}(0.0,2.0,0.0)
Vc4 = reverse(Vc2)
Vc4 .+= v_shift

V_curve_x0_2 = reverse(deepcopy(V_curve_x0))
V_curve_x0_2 .+= v_shift
V_curve_y0 = reverse(VT[ind_curve_y[2:end-1]])
V_curve_y0_2 = deepcopy(V_curve_y0)
v_shift = Point{3,Float64}(0.0,2.0,0.0)
V_curve_y0_2 .+= v_shift

# Meshing bottom surface 
Vm_region_curve = [Vc1; Vc2; Vc3; Vc4]
Fm,Vm,Cm = regiontrimesh((Vm_region_curve,),([1],),(pointSpacingSide))
invert_faces!(Fm)

Fb_tet,Vb_tet,Cb_tet  = joingeom(FT,VT,Ff,Vf,Fb,Vb,Fm,Vm)
Fb_tet,Vb_tet = mergevertices(Fb_tet,Vb_tet)

if element_formulation == :quadratic
    element_type = Tet10{Int}
elseif element_formulation == :linear 
    element_type = Tet4{Int}    
end

stringOpt = "paAqY"
vol1 = tetrahedral_volume_factor*(pointSpacing^3 / (6.0*sqrt(2.0)))
E_tet,V_tet,CE_tet,Fb_tet,Cb_tet = tetgenmesh(Fb_tet,Vb_tet; facetmarkerlist=Cb_tet, region_vol=vol1, stringOpt, element_type=element_type)

F_tet = element2faces(E_tet) # Triangular faces
CE_F_tet = repeat(CE_tet,inner=4)

## Join and merge pentahedral and tetrahedral 
V = [V_tet; V_penta]

Q = RotXYZ(-pi/2.0, 0.0, 0.0)
# Q2 = RotXYZ(0.0, -pi/2.0, 0.0)
# Q = Q2*Q1
V = [Q*Point{3, Float64}(v[1]-axialLength/2.0, v[2], v[3]) for v ∈ V] 

indexShift = length(V_tet)
E_penta = [e.+indexShift for e in E_penta]

V, indUnique, indMap = mergevertices(V; pointSpacing=pointSpacing)

indexmap!(E_penta, indMap)
indexmap!(E_tet, indMap)
indexmap!(F_tet, indMap)
indexmap!(Fb_tet, indMap)

FE_penta = element2faces(E_penta)

## Define Abaqus INP file
jobName = "Auxetic_structure"
partName_1 = "Auxetic_structure"
instanceName_1 = "Auxetic_structure-1"

sectionName_1 = "Substrate"
if element_formulation == :quadratic
    elementType_1 = "C3D10H" 
elseif element_formulation == :linear
    elementType_1 = "C3D4H"
end
nodeSetName_1 = "NodeSet-1_substrate"
elementSetName_1 = "ElementSet-1_substrate"
materialName_1 = "Mat_substrate"

sectionName_2 = "Film"
if element_formulation == :quadratic
    elementType_2 = "C3D15"# "C3D6"
elseif element_formulation == :linear
    elementType_2 = "C3D6"
end
nodeSetName_2 = "NodeSet-2_film"
elementSetName_2 = "ElementSet-2_film"
materialName_2 = "Mat_film"

nodeSetName_3 = "NodeSet-3_Bottom"
nodeSetName_4 = "NodeSet-4_Front"
nodeSetName_5 = "NodeSet-5_Back"

nodeSetName_assembly_REF = "REF_X"
nodeSetName_assembly_back = "NodeSet-6_MidBack"
nodeSetName_assembly_front = "NodeSet-6_MidFront"

FE_penta_quad = FE_penta[2]
indBoundary_penta_quads = boundaryfaceindices(FE_penta_quad)
FE_penta_quad_boundary = FE_penta_quad[indBoundary_penta_quads]

Xp = [v[2] for v in V]
indBottomQuads = Vector{Int}()
for (i,f) in enumerate(FE_penta_quad_boundary)
    if all(Xp[f] .< tolLevel)
        push!(indBottomQuads,i)
    end
end
indBottomNodesTri = reduce(vcat,Fb_tet[Cb_tet.==4])
indBottomNodesQuad = reduce(vcat,FE_penta_quad_boundary[indBottomQuads])
indBottomNodes = unique([indBottomNodesTri;indBottomNodesQuad])

# ---------------------------------------------------------------------------
# Periodic end-node selection and one-to-one geometric pairing
#
# IMPORTANT:
# `unique(...)` preserves the order in which nodes happen to be encountered in
# the face-connectivity arrays. The two end faces generally have different
# connectivity orderings, so independently enumerating their node lists does
# not create valid periodic pairs.
#
# Select both end planes geometrically, then match every x-minus node to the
# unused x-plus node having the same transverse (y,z) coordinates.
# ---------------------------------------------------------------------------

x0 = -axialLength/2.0
x1 =  axialLength/2.0
endPlaneTol = tolLevel
pairTol = max(1.0e-8, pointSpacing*1.0e-5)

indFrontNodes_raw = findall(v -> abs(v[1] - x0) <= endPlaneTol, V)
indBackNodes_raw  = findall(v -> abs(v[1] - x1) <= endPlaneTol, V)

function pair_periodic_nodes(V, nodes_x0, nodes_x1; pairTol)
    nodes_x0 = sort(unique(collect(nodes_x0)); by=i -> (V[i][2], V[i][3]))
    nodes_x1 = unique(collect(nodes_x1))

    length(nodes_x0) == length(nodes_x1) || error(
        "Periodic end faces contain different node counts: " *
        "x0=$(length(nodes_x0)), x1=$(length(nodes_x1))."
    )

    available = trues(length(nodes_x1))
    nodes_x1_paired = Vector{Int}(undef, length(nodes_x0))
    pairErrors = Vector{Float64}(undef, length(nodes_x0))

    for (m, i0) in enumerate(nodes_x0)
        y0 = V[i0][2]
        z0 = V[i0][3]

        best_j = 0
        best_d2 = Inf

        @inbounds for j in eachindex(nodes_x1)
            available[j] || continue

            i1 = nodes_x1[j]
            dy = V[i1][2] - y0
            dz = V[i1][3] - z0
            d2 = dy*dy + dz*dz

            if d2 < best_d2
                best_d2 = d2
                best_j = j
            end
        end

        best_j == 0 && error("No unused periodic partner found for node $i0.")

        d = sqrt(best_d2)
        d <= pairTol || error(
            "Periodic pairing failed for x0 node $i0. " *
            "Nearest transverse mismatch is $d, but pairTol=$pairTol."
        )

        nodes_x1_paired[m] = nodes_x1[best_j]
        pairErrors[m] = d
        available[best_j] = false
    end

    all(.!available) || error("Some x1 periodic nodes were not paired.")

    return nodes_x0, nodes_x1_paired, pairErrors
end

indFrontNodes, indBackNodes, periodicPairErrors =
    pair_periodic_nodes(V, indFrontNodes_raw, indBackNodes_raw; pairTol=pairTol)

axialPairErrors = [
    abs((V[indBackNodes[m]][1] - V[indFrontNodes[m]][1]) - axialLength)
    for m in eachindex(indFrontNodes)
]

maximum(axialPairErrors) <= endPlaneTol || error(
    "Periodic end-node pairs do not span the expected axial length."
)

println("Periodic node pairs: ", length(indFrontNodes))
println("Maximum transverse pairing error: ", maximum(periodicPairErrors))
println("Maximum axial pairing error: ", maximum(axialPairErrors))


_, indMidBackNode = findmin(norm.(V.-P_back))
_, indMidFrontNode = findmin(norm.(V.-P_front))

file_io = open(fileName_inp, "w")

addHeader(file_io,jobName)
addPart(file_io, partName_1; firstTime = true)
    addNodes(file_io, V)
    addElements(file_io,E_tet,elementType_1; indexOffset=0)
    addElements(file_io,E_penta,elementType_2; indexOffset=length(E_tet))
    addIndexSet(file_io, nodeSetName_1, unique(reduce(vcat,E_tet)); type=:nodes, indexOffset=0)
    addIndexSet(file_io, nodeSetName_2, unique(reduce(vcat,E_penta)); type=:nodes, indexOffset=0)
    addIndexSet(file_io, nodeSetName_3, indBottomNodes; type=:nodes, indexOffset=0)
    addIndexSet(file_io, nodeSetName_4, indFrontNodes; type=:nodes, indexOffset=0)
    addIndexSet(file_io, nodeSetName_5, indBackNodes; type=:nodes, indexOffset=0)    
    addIndexSet(file_io, elementSetName_1, 1:length(E_tet); type=:elements)
    addIndexSet(file_io, elementSetName_2, 1:length(E_penta); type=:elements, indexOffset=length(E_tet))
    addSolidSection(file_io, elementSetName_1, materialName_1)
    addSolidSection(file_io, elementSetName_2, materialName_2)
endPart(file_io)

startAssembly(file_io; name="Assembly-1")
    addInstance(file_io; name=instanceName_1, part=partName_1)
    addNodes(file_io, [P_ref]) # Add reference node
    addIndexSet(file_io, nodeSetName_assembly_REF, [1]; type=:nodes, indexOffset=0)

    addIndexSet(file_io, nodeSetName_assembly_back, indMidBackNode; type=:nodes, indexOffset=0, instance=instanceName_1)
    addIndexSet(file_io, nodeSetName_assembly_front, indMidFrontNode; type=:nodes, indexOffset=0, instance=instanceName_1)

    for (m,i) in enumerate(indFrontNodes)        
        nodeSetName_NOW = "NODE_X0_$m"
        addIndexSet(file_io, nodeSetName_NOW, [i]; type=:nodes, indexOffset=0, instance=instanceName_1)
    end

    for (m,i) in enumerate(indBackNodes)        
        nodeSetName_NOW = "NODE_X1_$m"
        addIndexSet(file_io, nodeSetName_NOW, [i]; type=:nodes, indexOffset=0, instance=instanceName_1)
    end

    for m in eachindex(indFrontNodes)        
        constraint(file_io; constraintname="Equation_" * "X$m" * "_D1")
        equation(file_io; n=3, node_sets=["NODE_X1_$m", "NODE_X0_$m", nodeSetName_assembly_REF], vals=[[1, 1], [1, -1], [1, -1]])        

        constraint(file_io; constraintname="Equation_" * "X$m" * "_D2")
        equation(file_io; n=3, node_sets=["NODE_X1_$m", "NODE_X0_$m", nodeSetName_assembly_REF], vals=[[2, 1], [2, -1], [2, -1]])

        constraint(file_io; constraintname="Equation_" * "X$m" * "_D3")
        equation(file_io; n=3, node_sets=["NODE_X1_$m", "NODE_X0_$m", nodeSetName_assembly_REF], vals=[[3, 1], [3, -1], [3, -1]])
    end
    
endAssembly(file_io)

addMaterial(file_io; name=materialName_2, category="Hyperelastic", parameters=[110., 0.1, 0.5], user="", type="COMPRESSIBLE", properties="3")

addMaterial(file_io; name=materialName_1, category="Hyperelastic, mooney-rivlin", parameters=[0.1, 0.4, 0.00200067])

addBoundary(file_io; setName=instanceName_1*"."*nodeSetName_3, flag="YSYMM")

addBoundary(file_io; setName=nodeSetName_assembly_back, op="NEW", vals=[1, 1])
addBoundary(file_io; setName=nodeSetName_assembly_back, vals=[3, 3], skipHeading=true)  

# addBoundary(file_io; setName=nodeSetName_assembly_back, op="NEW", load_case=2, vals=[1, 1])
# addBoundary(file_io; setName=nodeSetName_assembly_back, vals=[3, 3], skipHeading=true)  

startStep(file_io; name="Prestretch", nlgeom="YES", type="Static", parameters="0.001, 1., 1e-15, 1.", inc=1000)
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[1, 1, 0.145])
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[2, 2], skipHeading=true)
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[3, 3], skipHeading=true)
    S = [   "** OUTPUT REQUESTS", 
            "**",
            "*Restart, write, frequency=0",
            "**", 
            "** FIELD OUTPUT: F-Output-1",
            "**", 
            "*Output, field, variable=PRESELECT",
            "**", 
            "** HISTORY OUTPUT: H-Output-1",
            "**", 
            "*Output, history, variable=PRESELECT"]
    addFree(file_io, S)
endStep(file_io)

startStep(file_io; name="Buckle", nlgeom="NO", type="Buckle", parameters="1, , 2, 1000", perturbation="")

    addBoundary(file_io; setName=nodeSetName_assembly_back, op="NEW", load_case=1, vals=[1, 1])
    addBoundary(file_io; setName=nodeSetName_assembly_back, vals=[3, 3], skipHeading=true)  
  
    addBoundary(file_io; setName=nodeSetName_assembly_back, op="NEW", load_case=2, vals=[1, 1])
    addBoundary(file_io; setName=nodeSetName_assembly_back, vals=[3, 3], skipHeading=true)  

    addBoundary(file_io; setName=instanceName_1*"."*nodeSetName_3, flag="YSYMM", op="NEW", load_case=1)
    addBoundary(file_io; setName=instanceName_1*"."*nodeSetName_3, flag="YSYMM", op="NEW", load_case=2)
    
    addBoundary(file_io; setName=nodeSetName_assembly_REF, op="NEW", load_case=1, vals=[1, 1, 1.45])
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[2, 2], skipHeading=true)
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[3, 3], skipHeading=true)

    addBoundary(file_io; setName=nodeSetName_assembly_REF, op="NEW", load_case=2, vals=[1, 1, 1.45])
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[2, 2], skipHeading=true)
    addBoundary(file_io; setName=nodeSetName_assembly_REF, vals=[3, 3], skipHeading=true)

    S = [   "**",  
            "** OUTPUT REQUESTS",
            "**",
            "*Restart, write, frequency=0",
            "**", 
            "** FIELD OUTPUT: F-Output-2",
            "**", 
            "*Output, field, variable=PRESELECT"]
    addFree(file_io, S)
endStep(file_io)

close(file_io)

# # -----------------------------------------------------------------------------
# Visualization


F_unit_p,V_unit_p = separate_vertices(F_unit,V_unit) # Give each face its own point set 
C_unit_p = simplex2vertexdata(F_unit_p,C_unit) # Convert face color data to vertex color data 

FTp,VTp = separate_vertices(FT,VT) # Give each face its own point set 
CTp = simplex2vertexdata(FTp,CT) # Convert face color data to vertex color data 

Fb_tetp,Vb_tetp = separate_vertices(Fb_tet,V) # Give each face its own point set 
Cb_tetp = simplex2vertexdata(Fb_tetp,Cb_tet) # Convert face color data to vertex color data 

indB = boundaryfaceindices(F_tet)
F_tetp,V_tetp = separate_vertices(F_tet[indB],V)

fig1 = Figure(size = (1200,800), fontsize = 20)
ax1 = Axis(fig1[1, 1], aspect = DataAspect(), title="Periodic units")
# hp1 = meshplot!(ax1, F_unit_p, V_unit_p, strokewidth=strokewidth, color=C_unit_p, colormap=Makie.Categorical(Makie.Reverse(:Spectral)), shading=false)
# Colorbar(fig[1, 2], hp1)
hp1 = meshplot!(ax1, F_unit_p[C_unit.==5], V_unit_p, strokewidth=strokewidth1, color=:lightgray, shading=false)
hp2 = meshplot!(ax1, F_unit_p[C_unit.!=5], V_unit_p, strokewidth=strokewidth1, color=:steelblue1, shading=false)

ax21 = AxisGeom(fig1[1, 2][1, 1], title="Solid mesh", azimuth=0.82*pi, elevation=0.2*pi)
ax21.zoom_mult[]=0.7

# hp21 = meshplot!(ax21, Fb_tetp, Vb_tetp, strokewidth=strokewidth, color=Cb_tetp, strokecolor=:black, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp21 = meshplot!(ax21, Fb_tetp, Vb_tetp, strokewidth=strokewidth2, color=:lightgray, strokecolor=:black)
hp31 = meshplot!(ax21, FE_penta[1], V, strokewidth=strokewidth2, color=:steelblue1)    
hp41 = meshplot!(ax21, FE_penta[2], V, strokewidth=strokewidth2, color=:steelblue1)    

ax22 = AxisGeom(fig1[1, 2][2, 1], title="Solid mesh", azimuth=-(0.82*pi), elevation=-0.2*pi)
ax22.zoom_mult[]=0.7

# hp22 = meshplot!(ax22, Fb_tetp, Vb_tetp, strokewidth=strokewidth, color=Cb_tetp, strokecolor=:black, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))
hp22 = meshplot!(ax22, Fb_tetp, Vb_tetp, strokewidth=strokewidth2, color=:lightgray, strokecolor=:black)
hp32 = meshplot!(ax22, FE_penta[1], V, strokewidth=strokewidth2, color=:steelblue1)    
hp42 = meshplot!(ax22, FE_penta[2], V, strokewidth=strokewidth2, color=:steelblue1)    

Legend(fig1[1, 3], [hp22, hp32], ["Tetrahedral core", "Pentahedral layer"])

screen = display(GLMakie.Screen(), fig1)

# Colorbar(fig[1, 4], hp21)

fig2 = Figure(size = (1200,800))
ax31 = AxisGeom(fig2[1, 1][1, 1], title="Tri. mesh", azimuth=0.75*pi, elevation=0.1*pi)
ax31.zoom_mult[]=0.7
hp51 = meshplot!(ax31, Fb_tetp, Vb_tetp, strokewidth=strokewidth3, color=:lightgray)
hp61 = meshplot!(ax31, FE_penta[1], V, strokewidth=strokewidth3, color=:steelblue1)    
hp71 = meshplot!(ax31, FE_penta[2], V, strokewidth=strokewidth3, color=:steelblue1)    

hp81 = scatter!(ax31, V[indBottomNodes], markersize=markersize1, color = :red)
hp91 = scatter!(ax31, V[indFrontNodes], markersize=markersize1, color = :green)
hp101 = scatter!(ax31, V[indBackNodes], markersize=markersize1, color = :orange)

ax32 = AxisGeom(fig2[1, 1][1, 2], title="Tri. mesh", azimuth=-pi/4, elevation=-0.1*pi)
ax32.zoom_mult[]=0.7
hp51 = meshplot!(ax32, Fb_tetp, Vb_tetp, strokewidth=strokewidth3, color=:lightgray)
hp61 = meshplot!(ax32, FE_penta[1], V, strokewidth=strokewidth3, color=:steelblue1)    
hp71 = meshplot!(ax32, FE_penta[2], V, strokewidth=strokewidth3, color=:steelblue1)    

hp82 = scatter!(ax32, V[indBottomNodes], markersize=markersize1, color = :red)
hp92 = scatter!(ax32, V[indFrontNodes], markersize=markersize1, color = :green)
hp102 = scatter!(ax32, V[indBackNodes], markersize=markersize1, color = :orange)


scatter!(ax32, P_front, markersize=25, color = :purple)
scatter!(ax32, P_back, markersize=25, color = :purple)
scatter!(ax32, P_ref, markersize=25, color = :purple)

Legend(fig2[1, 2], [hp82, hp92, hp102], ["Bottom nodes", "Front nodes", "Back nodes"])

# hpm = mesh!(ax3,GeometryBasics.Mesh(V,Fb_tet[Cb_tet.==4]), color=:red,  transparency=true)
# scatter!(ax3,V[indBottomNodes],markersize=25, color = :red)

# # Setting up slicing
# VE  = simplexcenter(E_tet,V)
# ZE = [v[3] for v in VE]
# Z = [v[3] for v in V]
# zMax = maximum(Z)
# zMin = minimum(Z)
# numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F_tet,V)))

# stepRange = range(zMin,zMax,numSlicerSteps)
# hSlider = Slider(fig[2, :], range = stepRange, startvalue = mean(stepRange),linewidth=30)

# on(hSlider.value) do z 
#     B = ZE .<= z
#     indShow = findall(B)
#     if isempty(indShow)              
#         hp7.visible=false        
#     else        
#         hp7.visible=true
#         F_tet = element2faces(E_tet[indShow])        
#         indB = boundaryfaceindices(F_tet)                
#         F_tetp,V_tetp = separate_vertices(F_tet[indB],V)        
#         hp7[1] = GeometryBasics.Mesh(V_tetp,F_tetp)
#     end
# end
# hSlider.selected_index[]+=1
# slidercontrol(hSlider,ax3)

screen = display(GLMakie.Screen(), fig2)