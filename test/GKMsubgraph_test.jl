G = GKMproj_space(4)
S1 = GKMsubgraph_from_vertices(G, ["x_0", "x_1", "x_2"]) # valid GKMsubgraph
S2 = GKMsubgraph_from_edges(G, [Edge(1, 2), Edge(2, 3), Edge(3,4), Edge(4, 1)]) # also valid GKM subgraph, but incompatible with connection C (see below)

@test GKM_isValid(G)
@test GKM_isValidSubgraph(S1)
@test GKM_isValidSubgraph(S2)

C = build_GKM_connection(G)
@test GKM_isValidConnection(C)
@test isCompatible(S1, C)
@test !isCompatible(S2, C)