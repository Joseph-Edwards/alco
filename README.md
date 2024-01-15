# The ALCO Package for GAP

The **ALCO** package provides tools for algebraic combinatorics, most of which was written for **GAP** during the author's Ph.D. program. This package provides implementations in **GAP** of octonion algebras, Jordan algebras, and certain important integer subrings of those algebras. It also provides tools to compute the parameters of t-designs in spherical and projective spaces (modeled as manifolds of primitive idempotent elements in a simple Euclidean Jordan algebra). Finally, this package provides tools to explore octonion lattice constructions, including octonion Leech lattices. 

## Legal

The ALCO package provides tools for algebraic combinatorics in GAP. 

Copyright (C) 2024 Benjamin Nasmith

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Setup

1. Install [GAP](https://www.gap-system.org/Download/). The ALCO package was prepared using version 4.12. 
 
2. Clone this repository in your GAP installation as `c:/gap-4.XX.Y/pkg/alco`.

3. Open a GAP session and use the command `LoadPackage(“alco”);`  to import all the commands from this package.

## Example Session

The ALCO package allows users to construct the octonion arithmetic (integer ring). 
In the example below, we construct the octonion arithmetic and verify that the 
basis vectors define an E8 lattice relative to the inner product shown:

```
gap> LoadPackage("alco");
true
gap> A := OctonionArithmetic(Integers);
<algebra of dimension 8 over Integers>
gap> g := List(Basis(A), x -> List(Basis(A), y -> Norm(x+y) - Norm(x) - Norm(y)));;
gap> Display(g);
[ [   2,   0,  -1,   0,   0,   0,   0,   0 ],
  [   0,   2,   0,  -1,   0,   0,   0,   0 ],
  [  -1,   0,   2,  -1,   0,   0,   0,   0 ],
  [   0,  -1,  -1,   2,  -1,   0,   0,   0 ],
  [   0,   0,   0,  -1,   2,  -1,   0,   0 ],
  [   0,   0,   0,   0,  -1,   2,  -1,   0 ],
  [   0,   0,   0,   0,   0,  -1,   2,  -1 ],
  [   0,   0,   0,   0,   0,   0,  -1,   2 ] ]
gap> Determinant(g);
1
gap> IsGossetLatticeGramMatrix(g);
true
```
      
We can also construct simple Euclidean Jordan algebras, including the Albert
algebra:

```
gap> J := AlbertAlgebra(Rationals);
<algebra of dimension 27 over Rationals>
gap> SemiSimpleType(Derivations(Basis(J)));
"F4"
gap> AsList(Basis(J));
[ i1, i2, i3, i4, i5, i6, i7, i8, j1, j2, j3, j4, j5, j6, j7, j8, k1, k2, k3, k4, k5, k6,
  k7, k8, ei, ej, ek ]
gap> List(Basis(J), x -> Trace(x));
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 ]
gap> List(Basis(J), x -> Norm(x));
[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/2, 1/2, 1/2 ]
gap> List(Basis(J), x -> Determinant(x));
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> One(J);
ei+ej+ek
gap> Determinant(One(J));
1  
```

The ALCO package also provides tools to construct octonion lattices, including 
octonion Leech lattices.

```
gap> short_vectors := Set(ShortestVectors(GramMatrix(A),4).vectors, y -> LinearCombination(B
asis(A), y));;
gap> filt := Filtered(short_vectors, x -> x^2 + x + 2*One(x) = Zero(x));; 
gap> Length(filt);
576
gap> s := Random(filt);
a3+a4+a5+a7+a8
gap> gens := List(Basis(A), x -> x*[[s,s,0],[0,s,s],ComplexConjugate([s,s,s])]);; 
gap> gens := Concatenation(gens);;
gap> Length(gens);
24
gap> L := OctonionLatticeByGenerators(gens, One(A)*IdentityMat(3)/2);
<free left module over Integers, with 24 generators>
Time of last command: 433 ms
gap> IsLeechLatticeGramMatrix(GramMatrix(L));
true
```