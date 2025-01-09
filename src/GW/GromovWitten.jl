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

"""
Calculate the contribution of the decorated tree to the equivariant Gromov--Witten invariant.
**Warnings:**
1. This does not divide by |Aut(decorated tree)|.
2. This needs more testing with small examples like P^1 or P^2
3. This is is not optimized yet.
"""
function GWTreeContribution(
  dt::GW_decorated_tree,
  R::GKM_cohomology_ring,
  con::GKM_connection,
  classes::Vector{FreeModElem{QQMPolyRingElem}};
  check::Bool=true)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  
  if check
    @req length(dt.marks) == length(classes) "incompatible numbers of marked points and cohomology classes"
    @req R.gkm == dt.gkm "decorated tree and cohomolon_gy ring don't belong to the same GKM graph"
    @req con.gkm == dt.gkm "decorated tree and GKM connection don't belong to the same GKM graph"
  end

  C = R.coeffRing
  H = R.cohomRing
  gkm = R.gkm
  
  res = C(1)//C(1)

  for e in edges(dt.tree)
    de = dt.edgeMult[e]
    res = res * h(imageOf(e, dt), de, con, R) // de
  end

  for v in vertices(dt.tree)

    valv = degree(dt.tree, v)
    res = res * (eulerClass(imageOf(v, dt), R))^(valv - 1)

    tmpSum = C(0)//C(1)

    for v2 in all_neighbors(dt.tree, v)
      e = Edge(v,v2)
      wev = weightClass(imageOf(e, dt), R) // edgeMult(e, dt)
      res = res // wev
      tmpSum = tmpSum + 1//wev
    end

    res = res * tmpSum^( valv - 3 +  count(i -> i==v, dt.marks) )
  end

  for i in 1:length(dt.marks)
    v = imageOf(dt.marks[i], dt)
    pullback = classes[i][v] # pull back input class at i^th marked point to the fixed point v
    res = res * pullback
  end

  return res
end

"""
Calculate h(epsilon, d) as in [Liu--Sheshmani, Lemma 4.5, p. 16].
"""
function h(e::Edge, d::Int, con::GKM_connection, R::GKM_cohomology_ring; check::Bool=true)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  gkm = con.gkm
  C = R.coeffRing

  if check
    @req con.gkm == R.gkm "GKM connection and cohomology ring don't belong to the same GKM graph"
    @req has_edge(gkm.g, e) "edge not found in GKM graph"
    @req d>0 "d is non-positive"
  end

  we = weightClass(e, R) # weight of the edge e

  res = ( C((-1)^d) * ZZ(d)^(2d) ) // ( factorial(ZZ(d))^2 )
  res = res // ( (we)^(2d) )

  for v in all_neighbors(gkm.g, src(e))

    v == dst(e) && continue 

    ei = Edge(src(e), v)
    wei = weightClass(ei, R)
    ai = con.a[e][ei]
    res = res * b(1//d * we, wei, d*ai, C)
  end

  return res
end

"""
Calculate b(u,w,a) as in [Liu--Sheshmani, Lemma 4.5, p.16]. C is the coefficient ring
"""
function b(u::QQMPolyRingElem, w::QQMPolyRingElem, a::ZZRingElem, C::QQMPolyRing)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  res = C(1) // C(1) # make sure this has FracFieldElem type.
  if a >= 0
    for j in 0:a
      res = res // (w - j*u)
    end
  else 
    for j in 1:(-a-1)
      res = res // (w + j*u)
    end
  end
  return res
end