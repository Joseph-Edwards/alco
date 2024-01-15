LoadPackage("alco");

g := List(OctonionE8Basis, x -> List(OctonionE8Basis, y -> Norm(x+y) - Norm(x) - Norm(y)));
short := ShortestVectors(g, 4);
vecs := Set(short.vectors, x -> LinearCombination(OctonionE8Basis, x));
vecs := Union(vecs, - vecs);
filt := Filtered(vecs, x -> x^2 + x + 2*One(x) = Zero(x));
s := Random(filt);
sc := ComplexConjugate(s);

leech_basis := Concatenation(List(OctonionE8Basis, x -> x*[[sc,sc,0],[0,sc,sc],[s,s,s]]));

L := OctonionLatticeByGenerators(leech_basis, IdentityMat(3)*One(Oct)/2);

Display(IsLeechLatticeGramMatrix(GramMatrix(L)));

short := ShortestVectors(GramMatrix(L), 4);;
vecs := Set(short.vectors, x -> RealToOctonionVector(L, x));

commutative := Filtered(vecs, x -> IsCommutative(x));
projectors := Set(commutative, x -> VectorToIdempotentMatrix(x));
reflectors := Set(projectors, x -> One(x) - 2*x);
alb_refl := Set(reflectors, x -> HermitianMatrixToAlbertVector(x));

