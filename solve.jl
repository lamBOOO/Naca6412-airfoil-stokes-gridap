using Gridap
using GridapGmsh

mshfile = joinpath("NACA6412airfoil.msh")
model = GmshDiscreteModel(mshfile)

reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
reffeₚ = ReferenceFE(lagrangian,Float64,1)

labels = get_face_labeling(model)

V = TestFESpace(model,reffeᵤ,labels=labels,dirichlet_tags=["inflow", "airfoil", "wall"],conformity=:H1)
Q = TestFESpace(model,reffeₚ,conformity=:H1)
Y = MultiFieldFESpace([V,Q])

u_noslip = VectorValue(0,0)
u_inflow(x) = VectorValue([1.0-4.0*x[2]*x[2], x[1]*0.0])

U = TrialFESpace(V,[u_inflow,u_noslip,u_noslip])
P = TrialFESpace(Q)
X = MultiFieldFESpace([U,P])

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)

f = VectorValue(0.0,0.0)
a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ
l((v,q)) = ∫( v⋅f )dΩ

op = AffineFEOperator(a,l,X,Y)

uh, ph = solve(op)

writevtk(Ω,"results",order=2,cellfields=["uh"=>uh,"ph"=>ph]); @info "wrote results.vtk"
























plot = false

if plot
# PLOTTING PLAYGROUND (IGNORE IF YOU USE PARAVIEW)
# TOOD:
# - clean this up
using CairoMakie
using Gridap.CellData
using Gridap.Arrays


search_method = KDTreeSearch(num_nearest_vertices=5)
uhi = Interpolable(uh; searchmethod=search_method)
cache_u = return_cache(uhi, Gridap.Point(0.0, 0.0))
function uh_help(x)
  return evaluate!(cache_u, uhi, Gridap.Point(x))
end

function uh_ext(x)
  try
    return uh_help(x) |> x->[x[1], x[2]]
  catch
    return 1E-6*[1,1]
  end
end

uh_ext_pt(x) = Point2f(uh_ext(x))

begin
fig = Figure(
  # size = (800, 800)
)
ax = Axis(
  fig[1, 1]
  ,spinewidth=1.0
  ,title="stream lines"
  ,xlabel="x"
  ,ylabel="y"
  ,xminorticksvisible = true
  ,yminorticksvisible = true
  ,xticks = LinearTicks(6)
  ,yticks = LinearTicks(6)
  # , backgroundcolor = "black"
  , aspect = 1.25/0.5
  ,limits = (-1.25, 1.25, -0.5, 0.5)
)
streamplot!(ax,
  x->Point2f(uh_ext(x)), -1.25..1.25, -0.5..0.5, colormap = :grays
  # ,stepsizes = 0.1
  ,quality = 5
  ,arrow_size = 10
)
fig
# f = Figure()
# Axis(f[1, 1])
# xs = LinRange(0, 10, 100)
# ys = LinRange(0, 15, 100)
# zs = [cos(x) * sin(y) for x in xs, y in ys]
# contour!(zs,levels=-1:0.1:1)
# f

end


nodes = model.grid.node_coordinates

# xs = map(x->x[1], nodes)
# ys = map(x->x[2], nodes)

# uh_nodes = uh.(nodes)

# us = map(x->x[1], uh_nodes)
# vs = map(x->x[2], uh_nodes)
# # vs = [uh(VectorValue(x,y))[1] for x in xs, y in ys]
# strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
# # function method
# # arrow_fun(x) = Point2f(sin(x[1])*cos(x[2]), -cos(x[1])*sin(x[2]))
# arrows!(xs, ys, us, vs, arrowsize = 10, lengthscale = 0.3,
#     arrowcolor = strength, linecolor = strength)
# fig




# plot airfoild
airfoil_tags = model.face_labeling.tag_to_entities[findfirst(x->x=="airfoil",model.face_labeling.tag_to_name)]

tag_assignment = model.face_labeling.d_to_dface_to_entity[2]

all_faces = model.grid_topology.n_m_to_nface_to_mfaces[2]
airfoil_filter = [map(x -> tag_assignment[x[1]]==tag, enumerate(model.grid_topology.n_m_to_nface_to_mfaces[2])) for tag in airfoil_tags]

airfoil_faces_indices = map(x->[i[1] for i in enumerate(x) if i[2]], airfoil_filter)
airfoil_faces = [all_faces[i] for i in airfoil_faces_indices]

coords = model.grid.node_coordinates



airfoil_points = map(y->map(x->[coords[x[1]].data,coords[x[2]].data], y), airfoil_faces)
airfoil_points = map(x->unique(vcat(x...)), airfoil_points)

closed_pt_loops = [[] for i=1:length(airfoil_faces)]
for index in 1:length(airfoil_faces)
  closed_pt_loops[index] = [airfoil_faces[index][1][1]]
  for i=1:length(airfoil_faces[index])-1
    # @info closed_pt_loops[index]
    # @info "ho"
    nextpos = filter(x->x[1] == closed_pt_loops[index][end], airfoil_faces[index])
    # @show nextpos
    # nextneg = filter(x->x[2] == closed_pt_loops[index][end], airfoil_faces[index])
    # @show nextneg
    prevpos = filter(x->x[2] == closed_pt_loops[index][1], airfoil_faces[index])
    # @show prevpos
    # prevneg = filter(x->x[1] == closed_pt_loops[index][1], airfoil_faces[index])
    # @show prevneg
    # @info "hi"
    if length(nextpos) > 0
      push!(closed_pt_loops[index], nextpos[1][2])
    elseif length(prevpos) > 0
      pushfirst!(closed_pt_loops[index], prevpos[end][1])
    else
      error("no next prev")
    end
  end
end

airfoil_coords_ordered = [map(x->x.data, coords[closed_pt_loops[i]]) for i=1:length(airfoil_faces)]
airfoil_coords_combined = vcat(airfoil_coords_ordered...)
poly!(airfoil_coords_combined, color=:black)
fig

save("./result.pdf", fig)
save("./result.png", fig)
save("./result.svg", fig)


end
