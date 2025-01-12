"""
Create a decorated tree corresponding to the given GKM graph from the given data.
"""
function decoratedTree(
  gkm::AbstractGKM_graph,
  tree::Graph{Undirected},
  vDict::Dict{Int, Int},
  edgeMult::Dict{Edge, Int},
  marks::Vector{Int};
  check::Bool=true)::GW_decorated_tree

  nv = n_vertices(tree)

  for e in keys(edgeMult)
    edgeMult[reverse(e)] = edgeMult[e]
  end

  if check
    @req nv == length(vDict) "tree and vDict have different lengths"
    @req all([has_edge(gkm.g, Edge(vDict[src(e)], vDict[dst(e)])) for e in edges(tree)]) "image of tree edge does not exist in GKM graph"
    @req all([0 < edgeMult[e] for e in edges(tree)]) "non-positive edge multiplicity"
  end

  return GW_decorated_tree(gkm, tree, vDict, edgeMult, marks)
end

"""
Return the image of the edge e of the decorated tree t in the corresponding GKM graph.
"""
function imageOf(e::Edge, t::GW_decorated_tree)::Edge
  return Edge(t.vDict[src(e)], t.vDict[dst(e)])
end

"""
Return the image of the vertex v of the decorated tree t in the corresponding GKM graph
"""
function imageOf(v::Int, t::GW_decorated_tree)::Int
  return t.vDict[v]
end

function edgeMult(e::Edge, dt::GW_decorated_tree)::Int
  if e in keys(dt.edgeMult)
    return dt.edgeMult[e]
  elseif reverse(e) in keys(dt.edgeMult)
    return dt.edgeMukt[reverse(e)]
  end
  @req false "edge has no multiplicity assigned in decorated tree"
end
