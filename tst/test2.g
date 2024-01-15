LoadPackage("alco");
A := OctonionArithmetic(Integers);
g := List(Basis(A), x -> List(Basis(A), y -> Norm(x+y) - Norm(x) - Norm(y)));
short := ShortestVectors(g, 4);
vecs := Set(short.vectors, x -> LinearCombination(Basis(A), x));
vecs := Union(vecs, - vecs);
filt := Filtered(vecs, x -> x^2 + x + 2*One(x) = Zero(x));
s := Random(filt);
sc := ComplexConjugate(s);

frame := Filtered(filt, x -> x mod 2 = s mod 2);

leech_basis := Concatenation(List(Basis(A), x -> x*[[sc,sc,0],[0,sc,sc],[s,s,s]]));

S := f -> [
    [2,0,0],
    ComplexConjugate([f, 0, f]),
    [1, 1, f]
]*One(f); 

L := OctonionLatticeByGenerators(leech_basis, IdentityMat(3)*One(A)/2);

Display(IsLeechLatticeGramMatrix(GramMatrix(L)));

gens := Concatenation(S(s){[1,2]}, List(frame, f -> S(f)[3]));

reflectors := List(gens, x -> One(x)*IdentityMat(3)  - 2*VectorToIdempotentMatrix(x));

J := SimpleEuclideanJordanAlgebra(3, 8, Basis(A));

jord_refl := List(reflectors, x -> HermitianMatrixToJordanVector(x, J));

RandomClosure(jord_refl{[1..5]}, P, 5, 1);;
