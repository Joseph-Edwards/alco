#     This file belongs to ALCO: Tools for Algebraic Combinatorics.
#     Copyright (C) 2024, Benjamin Nasmith

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Quaternion tools

# THIS IS DIFFERENT THAT STANDARD GAP. Might need to remove. 
#   Perhaps I do not implement "Imaginary Part" for quaternions or octonions in my package. Encourage use of x - Re(x) in order to compute Im(x). 

# InstallMethod( RealPart,
#     "for a quaternion",
#     [ IsQuaternion and IsSCAlgebraObj ],
#     quat -> (1/2)*(quat + ComplexConjugate(quat)) 
#     );

# InstallMethod( ImaginaryPart,
#     "for a quaternion",
#     [ IsQuaternion and IsSCAlgebraObj ],
#     quat -> (1/2)*(quat - ComplexConjugate(quat)) 
#     );

InstallMethod( Trace, 
    "for a quaternion",
    [ IsQuaternion and IsSCAlgebraObj ],
    quat -> ExtRepOfObj(quat + ComplexConjugate(quat))[1]
    );

InstallMethod( Norm, 
    "for a quaternion",
    [ IsQuaternion and IsSCAlgebraObj ],
    quat -> ExtRepOfObj(quat*ComplexConjugate(quat))[1]
    );

InstallValue( QuaternionD4Basis, Basis(QuaternionAlgebra(Rationals), 
    List([  [  -1/2,  -1/2,  -1/2,   1/2 ],
            [  -1/2,  -1/2,   1/2,  -1/2 ],
            [  -1/2,   1/2,  -1/2,  -1/2 ],
            [     1,     0,     0,     0 ] ], x -> ObjByExtRep(FamilyObj(One(QuaternionAlgebra(Rationals))),x)
        )) );

# Icosian and Golden Field Tools

InstallValue( sigma, (1-Sqrt(5))/2 );

InstallValue( tau, (1+Sqrt(5))/2 );

# Previous versions of these functions:

# InstallGlobalFunction(GoldenRationalComponent, function(z)
#     return Trace(Field(sigma),Rationals,z)/2;
# end );

# InstallGlobalFunction(GoldenIrrationalComponent, function(z)
#     return (z-Trace(Field(sigma),Rationals,z)/2)/Sqrt(5 );
# end );

InstallGlobalFunction( GoldenModSigma, function(q)
    # Check that q belongs to the quadratic field containing Sqrt(5).
    if not q in NF(5,[1,4]) then return fail; fi;
    # Compute the coefficients in the basis [1,sigma] and return the 1 coefficient.
    return Coefficients(Basis(NF(5,[1,4]), [1, sigma]), q)[1];
end );

InstallValue( IcosianH4Basis, Basis(QuaternionAlgebra(Field(Sqrt(5))),
    List([  [ 0, -1, 0, 0 ], 
            [ 0, -1/2*E(5)^2-1/2*E(5)^3, 1/2, -1/2*E(5)-1/2*E(5)^4 ],
            [ 0, 0, -1, 0 ], 
            [ -1/2*E(5)-1/2*E(5)^4, 0, 1/2, -1/2*E(5)^2-1/2*E(5)^3 ] ], 
        x -> ObjByExtRep(FamilyObj(One(QuaternionAlgebra(Field(Sqrt(5))))),x)
        )) );

# Function to construct an octonion algebra in a standard basis (i.e. e[1]*e[2] = e[4], and cycle indices). 
InstallGlobalFunction( OctonionAlgebra, function(F)
    local T, alg;
    # Verify the input.
    if not (IsField(F) or IsIntegers(F)) then return fail; fi;
    # Precalculated structure constants table. 
    T :=  [ [ [ [ 8 ], [ -1 ] ], [ [ 4 ], [ 1 ] ], [ [ 7 ], [ 1 ] ], [ [ 2 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], 
                [ [ 5 ], [ -1 ] ], [ [ 3 ], [ -1 ] ], [ [ 1 ], [ 1 ] ] ], 
            [ [ [ 4 ], [ -1 ] ], [ [ 8 ], [ -1 ] ], [ [ 5 ], [ 1 ] ], [ [ 1 ], [ 1 ] ], [ [ 3 ], [ -1 ] ], 
                [ [ 7 ], [ 1 ] ], [ [ 6 ], [ -1 ] ], [ [ 2 ], [ 1 ] ] ], 
            [ [ [ 7 ], [ -1 ] ], [ [ 5 ], [ -1 ] ], [ [ 8 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], [ [ 2 ], [ 1 ] ], 
                [ [ 4 ], [ -1 ] ], [ [ 1 ], [ 1 ] ], [ [ 3 ], [ 1 ] ] ], 
            [ [ [ 2 ], [ 1 ] ], [ [ 1 ], [ -1 ] ], [ [ 6 ], [ -1 ] ], [ [ 8 ], [ -1 ] ], [ [ 7 ], [ 1 ] ], 
                [ [ 3 ], [ 1 ] ], [ [ 5 ], [ -1 ] ], [ [ 4 ], [ 1 ] ] ], 
            [ [ [ 6 ], [ -1 ] ], [ [ 3 ], [ 1 ] ], [ [ 2 ], [ -1 ] ], [ [ 7 ], [ -1 ] ], [ [ 8 ], [ -1 ] ], 
                [ [ 1 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [ 5 ], [ 1 ] ] ], 
            [ [ [ 5 ], [ 1 ] ], [ [ 7 ], [ -1 ] ], [ [ 4 ], [ 1 ] ], [ [ 3 ], [ -1 ] ], [ [ 1 ], [ -1 ] ], 
                [ [ 8 ], [ -1 ] ], [ [ 2 ], [ 1 ] ], [ [ 6 ], [ 1 ] ] ], 
            [ [ [ 3 ], [ 1 ] ], [ [ 6 ], [ 1 ] ], [ [ 1 ], [ -1 ] ], [ [ 5 ], [ 1 ] ], [ [ 4 ], [ -1 ] ], 
                [ [ 2 ], [ -1 ] ], [ [ 8 ], [ -1 ] ], [ [ 7 ], [ 1 ] ] ], 
            [ [ [ 1 ], [ 1 ] ], [ [ 2 ], [ 1 ] ], [ [ 3 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [ 5 ], [ 1 ] ], 
                [ [ 6 ], [ 1 ] ], [ [ 7 ], [ 1 ] ], [ [ 8 ], [ 1 ] ] ], 0, Zero(F) ];
    # Define the algebra and properties.
    alg := AlgebraByStructureConstantsArg([F, T, "e"], IsSCAlgebraObj and IsOctonion );
    SetFilterObj( alg, IsOctonionAlgebra );
    SetGramMatrix( alg, 2*IdentityMat(8) );
    return alg;
end );

InstallMethod( Trace,
    "for an octonion",
    [ IsOctonion and IsSCAlgebraObj ],
    oct -> ExtRepOfObj(oct)*GramMatrix(FamilyObj(oct)!.fullSCAlgebra)*ExtRepOfObj(One(oct)) 
    );

InstallMethod( Norm,
    "for an octonion",
    [ IsOctonion and IsSCAlgebraObj ],
    oct -> ExtRepOfObj(oct)*GramMatrix(FamilyObj(oct)!.fullSCAlgebra)*ExtRepOfObj(oct)/2 
    );

InstallMethod( Norm,
    "for an octonion list",
    [ IsOctonionCollection ],
    function(vec)
        return Sum(List(vec, x -> Norm(x)) );
    end );

InstallMethod( ComplexConjugate,
    "for an octonion",
    [ IsOctonion and IsSCAlgebraObj ],
    oct -> One(oct)*Trace(oct) - oct
    );  

InstallMethod( RealPart,
    "for an octonion",
    [ IsOctonion and IsSCAlgebraObj ],
    oct -> (1/2)*Trace(oct)*One(oct)
    );

# Avoiding use of ImaginaryPart for octonions since the built in GAP quaternion version of this command is easily misunderstood. 

# InstallMethod( ImaginaryPart,
#     "for an octonion",
#     [ IsOctonion and IsSCAlgebraObj ],
#     oct -> oct - RealPart(oct)
#     );

InstallValue( Oct, OctonionAlgebra(Rationals) );  

# Octonion arithmetic tools

InstallValue( OctonionE8Basis, Basis(Oct,
    List(
        [[ -1/2, 0, 0, 0, 1/2, 1/2, 1/2, 0 ], 
    [ -1/2, -1/2, 0, -1/2, 0, 0, -1/2, 0 ], 
    [ 0, 1/2, 1/2, 0, -1/2, 0, -1/2, 0 ], 
    [ 1/2, 0, -1/2, 1/2, 1/2, 0, 0, 0 ],
    [ 0, -1/2, 1/2, 0, -1/2, 0, 1/2, 0 ], 
    [ 0, 1/2, 0, -1/2, 1/2, -1/2, 0, 0 ], 
    [ -1/2, 0, -1/2, 1/2, -1/2, 0, 0, 0 ], 
    [ 1/2, 0, 0, -1/2, 0, 1/2, 0, -1/2 ] ], 
        x -> ObjByExtRep(FamilyObj(One(Oct)),x)
        )) );       

InstallGlobalFunction( OctonionArithmetic, function(F, option...)
    local T, alg, basis;
    if not (IsField(F) or F = Integers) then return fail; fi;
    if Length(option) = 0 then 
        basis := OctonionE8Basis;
    else 
        basis := Basis(UnderlyingLeftModule(OctonionE8Basis), OctonionE8Basis*Inverse(OctonionE8Basis[8]) );
    fi;
    T := StructureConstantsTable(basis );
    alg := AlgebraByStructureConstantsArg([F, T, "a"], IsSCAlgebraObj and IsOctonionArithmeticElement );
    SetFilterObj( alg, IsOctonionAlgebra );
    SetGramMatrix( alg, 
      [ [   2,   0,  -1,   0,   0,   0,   0,   0 ],
        [   0,   2,   0,  -1,   0,   0,   0,   0 ],
        [  -1,   0,   2,  -1,   0,   0,   0,   0 ],
        [   0,  -1,  -1,   2,  -1,   0,   0,   0 ],
        [   0,   0,   0,  -1,   2,  -1,   0,   0 ],
        [   0,   0,   0,   0,  -1,   2,  -1,   0 ],
        [   0,   0,   0,   0,   0,  -1,   2,  -1 ],
        [   0,   0,   0,   0,   0,   0,  -1,   2 ] ] 
    );      
    return alg;
end );

InstallMethod( \mod, 
    "For an octonion integer", 
    [IsOctonionArithmeticElement, IsPosInt], 0,
    function(a,m)
        return ObjByExtRep(FamilyObj(One(a)), List(ExtRepOfObj(a), i -> i mod m) );
    end );

InstallGlobalFunction( OctonionToRealVector,
    function(basis, x) 
    # In the case of an octonion basis.
    if IsBasis(basis) and IsOctonionCollection(UnderlyingLeftModule(basis)) and IsOctonionCollection(x) then 
        return Flat(List(x, y -> Coefficients(basis, y)) );
    fi;
    # In the case of an octonion lattice.
    if IsOctonionLattice(basis) and IsOctonionCollection(x) then 
        return SolutionMat(LLLReducedBasisCoefficients(basis),
            OctonionToRealVector(CanonicalBasis(UnderlyingOctonionRing(basis)), x)
        );
    fi;
    return fail;
    end );

InstallGlobalFunction( RealToOctonionVector,
    function(basis, x)
        local n, temp;
        if Length(x) mod 8 <> 0 or not IsHomogeneousList(x) then 
            return fail; 
        fi;
        # In the case of an octonion basis,
        if IsBasis(basis) and IsOctonionCollection(UnderlyingLeftModule(basis)) then
            n := Length(x)/8;
            temp := List([1..n], m -> x{[(m-1)*8+1 .. m*8]} );;
            return List(temp, y -> LinearCombination(basis, y) );
        fi;
        # In the case of an octonion lattice.
        if IsOctonionLattice(basis) then 
            temp := LinearCombination(LLLReducedBasisCoefficients(basis), x );
            temp := RealToOctonionVector(CanonicalBasis(UnderlyingOctonionRing(basis)), temp );
            return temp;
        fi;
        return fail;
    end );

InstallGlobalFunction( VectorToIdempotentMatrix, 
    function(x)
        local temp;
        if not IsHomogeneousList(x) or not IsAssociative(x) or not (IsCyc(x[1]) or IsQuaternion(x[1]) or IsOctonion(x[1])) then 
            return fail;
        fi;
        if IsAssociative(x) then 
            temp := TransposedMat([ComplexConjugate(x)])*[x];
            return temp/Trace(temp );
        fi; 
        return fail;
    end );

InstallGlobalFunction( WeylReflection, 
    function(r,x)
        local R;
        R := VectorToIdempotentMatrix(r );
        if R = fail or not IsHomogeneousList(Flat([r,x])) then 
            return fail; 
        fi;
        return x - 2*x*R;
    end );


# Jordan Algebra Tools

InstallMethod( Rank, 
    "for a Jordan Algebra",
    [ IsJordanAlgebra ],
    J -> JordanRank(J)
    );

InstallMethod( Rank, 
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj ],
    j -> Rank(FamilyObj(j)!.fullSCAlgebra)
    );

InstallMethod( Degree, 
    "for a Jordan Algebra",
    [ IsJordanAlgebra ],
    J -> JordanDegree(J)
    );

InstallMethod( Degree, 
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj ],
    j -> Degree(FamilyObj(j)!.fullSCAlgebra)
    );

InstallMethod( Trace,
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj and IsSCAlgebraObj ],
    # j -> (Rank(j)/Dimension(FamilyObj(j)!.fullSCAlgebra))*Trace(AdjointMatrix(Basis(FamilyObj(j)!.fullSCAlgebra), j))
    j -> ExtRepOfObj(j)*JordanBasisTraces(FamilyObj(j)!.fullSCAlgebra) 
    );

InstallMethod( Norm,
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj and IsSCAlgebraObj ],
    j -> Trace(j^2)/2
    );

InstallGlobalFunction( GenericMinimalPolynomial, 
    function(x)
        local p, prx, p1xr, r, j, qjx;
        p := [-Trace(x)];
        p1xr := [-Trace(x)];
        r := 1;
        repeat 
            r := r + 1;
            Append(p1xr, [-Trace(x^r)] ); 
            prx := (1/r)*(p1xr[r] + Sum(List([1..r-1], j -> 
                p1xr[r-j]*p[j]
            )) );
            Append(p, [prx] );
        until r = Rank(x );
        p := Reversed(p );
        Append(p, [1] );
        return p;
    end );

InstallMethod( Determinant,
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj and IsSCAlgebraObj ],
    j -> ValuePol(GenericMinimalPolynomial(j), 0)*(-1)^Rank(j) 
    );

InstallMethod( JordanAdjugate, 
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj ],
    function(j) return (-1)^(1+Rank(j))*ValuePol(ShiftedCoeffs(GenericMinimalPolynomial(j), -1), j ); end
    );

InstallMethod( IsPositiveDefinite, 
    "for a Jordan algebra element",
    [ IsJordanAlgebraObj ],
    function(j) local temp; 
        temp := GenericMinimalPolynomial(j ); 
        if 0 in temp then return false; fi;
        temp := Reversed(temp );
        temp := List([0..Rank(j)], n -> temp[n+1]*(-1)^n > 0 );
        if Set(temp) = [true] then return true; fi;
        return false;
    end );

# Function to construct a Hermitian simple Euclidean Jordan algebra basis.
InstallGlobalFunction( HermitianJordanAlgebraBasis , function(rho, comp_alg_basis)
    local peirce, d, F, conj, mat, frame, realbasis;
    # Ensure that the rank and degree are correct.
    if not (IsInt(rho) and rho > 1 and IsBasis(comp_alg_basis)) then return fail; fi;
    # Ensure that the second argument is the basis for a composition algebra. 
    d := Dimension(UnderlyingLeftModule(comp_alg_basis) );
    if not IsBasis(comp_alg_basis) and d in [1,2,4,8] then return fail; fi;
    # Ensure that for d = 8 the rank is appropriate.
    if d = 8 and rho > 3 then return fail; fi;
    # If not a quaternion or octonion collection, then it must belong to a field extensions or the rationals.
    if not IsOctonionCollection(comp_alg_basis) and not IsQuaternionCollection(comp_alg_basis) then 
        if LeftActingDomain(UnderlyingLeftModule(comp_alg_basis)) <> Rationals then return fail; fi;
        # Ensure the quadratic extension case is imaginary.
        if d = 2 and Set(comp_alg_basis, x -> RealPart(x) = x ) = [true] then return fail; fi;
    fi;
    # Define a function to construct a single entry Hermitian matrix.
    mat := function(n,x)
        local temp; 
        if Length(n) <> 2 then return fail; fi;
        temp := Zero(x)*IdentityMat(rho );
        temp[n[1]][n[2]] := x;
        temp[n[2]][n[1]] := ComplexConjugate(x );
        return temp;
    end;
    # Construct the diagonal matrices. 
    frame := Concatenation(List(IdentityMat(rho), x -> List([One(UnderlyingLeftModule(comp_alg_basis))], r -> DiagonalMat(x)*r)) );
    # Construct the off-diagonal matrices.
    peirce := List(Combinations([1..rho],2), n -> List(comp_alg_basis, x -> mat(n,x)) );
    # Return the matrices. 
    return Concatenation(frame, Concatenation(peirce) );
end );

# Function to convert a Hermitian matrix to a vector, or coefficients of a vector, in a Jordan algebra. 
InstallGlobalFunction( HermitianMatrixToJordanCoefficients, function(mat, comp_alg_basis)
    local temp, result, i, z, basis, realbasis, pos;
    # Verify that the input is a Hermitian matrix.
    if not IsMatrix(mat) and mat = TransposedMat(ComplexConjugate(mat)) then return fail; fi;
    # A record to record the results.
    result := rec( );
    # Ensure that the second input is a basis. 
    if not IsBasis(comp_alg_basis) then return fail; fi;
    # Record basic parameters. 
    basis := comp_alg_basis;
    result.d := Length(comp_alg_basis );
    result.rho := Length(DiagonalOfMat(mat) ); 
    # Define a basis for the diagonal components. 
    realbasis := Filtered(basis, x -> RealPart(x) = x );
    pos := List(realbasis, x -> Position(basis, x) );
    if Length(realbasis) = 0 then realbasis := [One(basis)]; fi;
    # Define a temporary list of coefficients.
    temp := []; 
    # Determine the coefficients due to the diagonal components of the matrix.
    # First find the canonical basis for the subalgebra spanned by the identity. 
    realbasis := Basis(Subalgebra(UnderlyingLeftModule(comp_alg_basis), [One(UnderlyingLeftModule(comp_alg_basis))]) );
    # Then find the coefficients.
    for i in [1..result.rho] do 
        Append(temp, Coefficients(realbasis, mat[i][i]) );
    od;
    # Find the off-diagonal coefficients next.
    for i in Combinations([1..result.rho],2) do 
        Append(temp, Coefficients(basis, mat[i[1]][i[2]]) );
    od;
    # Return the coefficients.
    return temp;
end );

# Function to convert a Hermitian matrix into a Jordan algebra vector. 
InstallGlobalFunction( HermitianMatrixToJordanVector, function(mat, J)
    local temp;
    # Verify that the second argument is a Jordan algebra with an off-diagonal basis defined. 
    if not (IsJordanAlgebra(J) and HasJordanOffDiagonalBasis(J)) then return fail; fi;
    # Verify that the matrix entries belong to the algebra spanned by the off-diagonal basis. 
    if not IsSubset(UnderlyingLeftModule(JordanOffDiagonalBasis(J)), Flat(mat)) then return fail; fi;\
    # Compute the coefficients, if possible.
    temp := HermitianMatrixToJordanCoefficients(mat, JordanOffDiagonalBasis(J) );
    if temp = fail then return temp; fi;
    # Return the coefficients in J-vector form.
    return LinearCombination(Basis(J), temp );
end );

# Function to construct Jordan algebra of Hermitian type.
InstallGlobalFunction( HermitianSimpleJordanAlgebra, function(rho, comp_alg_basis, F...)
    local jordan_basis, T, temp, coeffs, K, algebra, filter, n, m, z, l;
    # Ensure inputs are correct.
    if not (IsInt(rho) and rho > 1 and IsBasis(comp_alg_basis)) then return fail; fi;
    if Length(F) > 1 then return fail; fi;
    # Ensure that the optional field argument contains the left acting domain of the basis. 
    K := Rationals;
    if Length(F) = 1 then 
        if not IsField(F[1]) or not Set(Basis(LeftActingDomain(UnderlyingLeftModule(comp_alg_basis))), y -> y in F[1]) then 
            return fail; 
        else 
            K := F[1];
        fi;
    fi;
    # Evaluate basis vectors:
    jordan_basis := HermitianJordanAlgebraBasis(rho, comp_alg_basis );
    if jordan_basis = fail then return jordan_basis; fi;
    # Define an empty structure constants table.
    T := EmptySCTable(Size(jordan_basis), Zero(LeftActingDomain(UnderlyingLeftModule(comp_alg_basis))) );
    # Compute the structure constants.
    for n in [1..Size(jordan_basis)] do 
        for m in [1..Size(jordan_basis)] do 
            z := (1/2)*(jordan_basis[n]*jordan_basis[m] + jordan_basis[m]*jordan_basis[n] );
            # z := (jordan_basis[n]*jordan_basis[m] + jordan_basis[m]*jordan_basis[n] );
            temp := HermitianMatrixToJordanCoefficients(z, comp_alg_basis );
            coeffs := [];
            for l in [1..Size(jordan_basis)] do 
                if temp[l] <> 0 then 
                    Append(coeffs, [temp[l], l] );
                fi;
            od;
            SetEntrySCTable( T, n, m, coeffs );
        od;
    od;
    # Construct the algebra.
    filter:= IsSCAlgebraObj and IsJordanAlgebraObj;
    algebra := AlgebraByStructureConstantsArg([K, T], filter );
    SetFilterObj( algebra, IsJordanAlgebra );
    # Assign various attributes to the algebra.
    SetJordanRank( algebra, rho );
    SetJordanDegree( algebra, Length(comp_alg_basis) );
    SetJordanMatrixBasis( algebra, jordan_basis );
    SetJordanOffDiagonalBasis( algebra, comp_alg_basis );
    SetJordanHomotopeVector( algebra, One(algebra) );
    SetJordanBasisTraces( algebra, List(Basis(algebra), 
        j -> (Rank(j)/Dimension(FamilyObj(j)!.fullSCAlgebra))*Trace(AdjointMatrix(Basis(FamilyObj(j)!.fullSCAlgebra), j)))
        );
    return algebra;
end );





InstallGlobalFunction( JordanSpinFactor,  function(gram_mat)
    local result, T, n, m, z, temp, coeffs, filter;
    if not IsMatrix(gram_mat) or Inverse(gram_mat) = fail then return fail; fi;
    result := rec( );
    result.F := Field(Flat(gram_mat) );
    result.rho := 2;
    result.d := DimensionsMat(gram_mat)[1]-1;
    # Construct the algebra.
    T := EmptySCTable(result.d + 2, Zero(0) );
    SetEntrySCTable(T, 1, 1, [1, 1] );
    for m in [2..result.d + 2] do 
        SetEntrySCTable(T, m, 1, [1, m] );
        SetEntrySCTable(T, 1, m, [1, m] );
    od;
    for n in [2..result.d + 2] do 
        for m in [2..result.d + 2] do
            SetEntrySCTable( T, n, m, [gram_mat[n-1][m-1], 1] );
        od;
    od;
    filter:= IsSCAlgebraObj and IsJordanAlgebraObj;
    result.algebra := AlgebraByStructureConstantsArg([result.F, T], filter );
    SetFilterObj( result.algebra, IsJordanAlgebra );
    SetJordanRank( result.algebra, result.rho );
    SetJordanDegree( result.algebra, result.d );
    SetJordanHomotopeVector( result.algebra, One(result.algebra) );
    SetJordanBasisTraces( result.algebra, List(Basis(result.algebra), 
        j -> (Rank(j)/Dimension(FamilyObj(j)!.fullSCAlgebra))*Trace(AdjointMatrix(Basis(FamilyObj(j)!.fullSCAlgebra), j)))
        );
    return result.algebra;
end );



InstallMethod( JordanAlgebraGramMatrix, 
    "for a Jordan algebra",
    [ IsJordanAlgebra ],
    j -> List(Basis(j), x -> List(Basis(j), y -> Trace(x*y)))
    );

InstallGlobalFunction( JordanHomotope , function(ring, u, label...)
    local result, temp, filter, T, n, m, z, l, coeffs;
    if not IsJordanAlgebra(ring) then return fail; fi;
    result := rec( );
    result.rho := Rank(ring );
    result.d := Degree(ring );
    result.F := LeftActingDomain(ring );
    T := EmptySCTable(Dimension(ring), Zero(result.F) );
    for n in Basis(ring) do 
        for m in Basis(ring) do 
            z := n*(u*m) + (n*u)*m - u*(n*m );
            temp := Coefficients(Basis(ring), z );
            coeffs := [];
            for l in [1..Dimension(ring)] do 
                if temp[l] <> Zero(temp[l]) then 
                    Append(coeffs, [temp[l], l] );
                fi;
            od;
            SetEntrySCTable( T, Position(Basis(ring), n), Position(Basis(ring), m), coeffs );
        od;
    od;
    filter:= IsSCAlgebraObj and IsJordanAlgebraObj;
    if Length(label) > 0 and IsString(label[1]) then 
        result.algebra := AlgebraByStructureConstantsArg([result.F, T, label[1]], filter );
    else 
        result.algebra := AlgebraByStructureConstantsArg([result.F, T], filter );
    fi;
    SetFilterObj( result.algebra, IsJordanAlgebra );
    SetJordanRank( result.algebra, result.rho );
    SetJordanDegree( result.algebra, result.d );
    if HasJordanMatrixBasis( result.algebra) then 
        SetJordanMatrixBasis( result.algebra, JordanMatrixBasis(ring) );
    fi;
    if HasJordanOffDiagonalBasis(result.algebra) then 
        SetJordanOffDiagonalBasis( result.algebra, JordanOffDiagonalBasis(ring) );
    fi;
    SetJordanHomotopeVector( result.algebra, u );
    SetJordanBasisTraces( result.algebra, List(Basis(result.algebra), 
        j -> (Rank(j)/Dimension(FamilyObj(j)!.fullSCAlgebra))*Trace(AdjointMatrix(Basis(FamilyObj(j)!.fullSCAlgebra), j)))
        );
    return result.algebra;
end );

InstallGlobalFunction( SimpleEuclideanJordanAlgebra, function(rho, d, args...)
    local temp;
    temp := rec( );
    if not (IsInt(d) and IsInt(rho)) then return fail; fi;
    if rho < 2 then return fail; fi; 
    if rho > 2 and not d in [1,2,4,8] then return fail; fi;
    if d = 8 and rho > 3 then return fail; fi;
    if rho = 2 then 
        if Length(args) = 0 then 
            return JordanSpinFactor(IdentityMat(d+1) );
        elif IsMatrix(args[1]) and DimensionsMat(args[1]) = [d+1, d+1] and Inverse(args[1]) <> fail then 
            return JordanSpinFactor(args[1] );
        elif d in [1,2,4,8] and IsBasis(args[1]) then 
            return HermitianSimpleJordanAlgebra(rho, args[1] );
        else 
            return fail;
        fi;
    fi;

    if Length(args) = 0 then
        if d = 8 then 
            return HermitianSimpleJordanAlgebra(rho, Basis(OctonionAlgebra(Rationals)) );
        elif d = 4 then 
            return HermitianSimpleJordanAlgebra(rho, Basis(QuaternionAlgebra(Rationals)) );
        elif d = 2 then 
            return HermitianSimpleJordanAlgebra(rho, Basis(CF(4), [1, E(4)]) );
        elif d = 1 then 
            return HermitianSimpleJordanAlgebra(rho, Basis(Rationals, [1]) );
        fi;
    elif IsBasis(args[1]) then 
        return HermitianSimpleJordanAlgebra(rho, args[1] );
    else 
        return fail;
    fi;
end );

# Albert Algebra tools

# Albert algebra by structure constants, computed elsewhere. 
InstallGlobalFunction(AlbertAlgebra, function(F)
    local T, alg, jordan_basis, i, j, k, e, temp;
    T := [ [ [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 24 ], [ -1/2 ] ], [ [ 20 ], [ -1/2 ] ], 
      [ [ 23 ], [ -1/2 ] ], [ [ 18 ], [ 1/2 ] ], [ [ 22 ], [ -1/2 ] ], [ [ 21 ], [ 1/2 ] ], 
      [ [ 19 ], [ 1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 12 ], [ 1/2 ] ], 
      [ [ 15 ], [ 1/2 ] ], [ [ 10 ], [ -1/2 ] ], [ [ 14 ], [ 1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 11 ], [ -1/2 ] ], [ [ 9 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 1 ], [ 1/2 ] ], 
      [ [ 1 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 20 ], [ 1/2 ] ], [ [ 24 ], [ -1/2 ] ], 
      [ [ 21 ], [ -1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 19 ], [ 1/2 ] ], [ [ 23 ], [ -1/2 ] ], 
      [ [ 22 ], [ 1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], [ [ 16 ], [ -1/2 ] ], 
      [ [ 13 ], [ 1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 15 ], [ 1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 10 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 2 ], [ 1/2 ] ], 
      [ [ 2 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 23 ], [ 1/2 ] ], [ [ 21 ], [ 1/2 ] ], 
      [ [ 24 ], [ -1/2 ] ], [ [ 22 ], [ -1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 20 ], [ 1/2 ] ], 
      [ [ 17 ], [ -1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 16 ], [ -1/2 ] ], [ [ 14 ], [ 1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 12 ], [ -1/2 ] ], 
      [ [ 9 ], [ 1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 3 ], [ 1/2 ] ], 
      [ [ 3 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 18 ], [ -1/2 ] ], [ [ 17 ], [ 1/2 ] ], 
      [ [ 22 ], [ 1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], [ [ 19 ], [ -1/2 ] ], 
      [ [ 21 ], [ 1/2 ] ], [ [ 20 ], [ -1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 9 ], [ -1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 15 ], [ 1/2 ] ], [ [ 11 ], [ 1/2 ] ], 
      [ [ 13 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 4 ], [ 1/2 ] ] ], [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 22 ], [ 1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 18 ], [ 1/2 ] ], [ [ 23 ], [ 1/2 ] ], 
      [ [ 24 ], [ -1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 20 ], [ -1/2 ] ], [ [ 21 ], [ -1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 11 ], [ 1/2 ] ], [ [ 10 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], 
      [ [ 16 ], [ -1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [ 12 ], [ 1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [ 5 ], [ 1/2 ] ], [ [ 5 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 21 ], [ -1/2 ] ], 
      [ [ 23 ], [ 1/2 ] ], [ [ 20 ], [ -1/2 ] ], [ [ 19 ], [ 1/2 ] ], [ [ 17 ], [ 1/2 ] ], 
      [ [ 24 ], [ -1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 22 ], [ -1/2 ] ], [ [ 13 ], [ 1/2 ] ], 
      [ [ 15 ], [ -1/2 ] ], [ [ 12 ], [ 1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 9 ], [ -1/2 ] ], 
      [ [ 16 ], [ -1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 14 ], [ -1/2 ] ], [ [  ], [  ] ], 
      [ [ 6 ], [ 1/2 ] ], [ [ 6 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [ 26, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [ 19 ], [ -1/2 ] ], 
      [ [ 22 ], [ -1/2 ] ], [ [ 17 ], [ 1/2 ] ], [ [ 21 ], [ -1/2 ] ], [ [ 20 ], [ 1/2 ] ], 
      [ [ 18 ], [ 1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], [ [ 11 ], [ 1/2 ] ], 
      [ [ 14 ], [ 1/2 ] ], [ [ 9 ], [ -1/2 ] ], [ [ 13 ], [ 1/2 ] ], [ [ 12 ], [ -1/2 ] ], 
      [ [ 10 ], [ -1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [  ], [  ] ], 
      [ [ 7 ], [ 1/2 ] ], [ [ 7 ], [ 1/2 ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [ 26, 27 ], [ 1, 1 ] ], [ [ 17 ], [ -1/2 ] ], 
      [ [ 18 ], [ -1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 20 ], [ -1/2 ] ], [ [ 21 ], [ -1/2 ] ], 
      [ [ 22 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], [ [ 24 ], [ 1/2 ] ], [ [ 9 ], [ -1/2 ] ], 
      [ [ 10 ], [ -1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 16 ], [ 1/2 ] ], [ [  ], [  ] ], 
      [ [ 8 ], [ 1/2 ] ], [ [ 8 ], [ 1/2 ] ] ], 
  [ [ [ 24 ], [ -1/2 ] ], [ [ 20 ], [ 1/2 ] ], [ [ 23 ], [ 1/2 ] ], [ [ 18 ], [ -1/2 ] ], 
      [ [ 22 ], [ 1/2 ] ], [ [ 21 ], [ -1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 17 ], [ -1/2 ] ], 
      [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 8 ], [ -1/2 ] ], [ [ 4 ], [ -1/2 ] ], 
      [ [ 7 ], [ -1/2 ] ], [ [ 2 ], [ 1/2 ] ], [ [ 6 ], [ -1/2 ] ], [ [ 5 ], [ 1/2 ] ], 
      [ [ 3 ], [ 1/2 ] ], [ [ 1 ], [ -1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [  ], [  ] ], 
      [ [ 9 ], [ 1/2 ] ] ], [ [ [ 20 ], [ -1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 21 ], [ 1/2 ] ], 
      [ [ 17 ], [ 1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 23 ], [ 1/2 ] ], [ [ 22 ], [ -1/2 ] ], 
      [ [ 18 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 4 ], [ 1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 5 ], [ -1/2 ] ], [ [ 1 ], [ -1/2 ] ], 
      [ [ 3 ], [ 1/2 ] ], [ [ 7 ], [ -1/2 ] ], [ [ 6 ], [ 1/2 ] ], [ [ 2 ], [ -1/2 ] ], 
      [ [ 10 ], [ 1/2 ] ], [ [  ], [  ] ], [ [ 10 ], [ 1/2 ] ] ], 
  [ [ [ 23 ], [ -1/2 ] ], [ [ 21 ], [ -1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 22 ], [ 1/2 ] ], 
      [ [ 18 ], [ 1/2 ] ], [ [ 20 ], [ -1/2 ] ], [ [ 17 ], [ 1/2 ] ], [ [ 19 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 7 ], [ 1/2 ] ], [ [ 5 ], [ 1/2 ] ], 
      [ [ 8 ], [ -1/2 ] ], [ [ 6 ], [ -1/2 ] ], [ [ 2 ], [ -1/2 ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 1 ], [ -1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 11 ], [ 1/2 ] ], [ [  ], [  ] ], 
      [ [ 11 ], [ 1/2 ] ] ], [ [ [ 18 ], [ 1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 22 ], [ -1/2 ] ], 
      [ [ 24 ], [ -1/2 ] ], [ [ 23 ], [ 1/2 ] ], [ [ 19 ], [ 1/2 ] ], [ [ 21 ], [ -1/2 ] ], 
      [ [ 20 ], [ -1/2 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 2 ], [ -1/2 ] ], [ [ 1 ], [ 1/2 ] ], [ [ 6 ], [ 1/2 ] ], [ [ 8 ], [ -1/2 ] ], 
      [ [ 7 ], [ -1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 5 ], [ 1/2 ] ], [ [ 4 ], [ -1/2 ] ], 
      [ [ 12 ], [ 1/2 ] ], [ [  ], [  ] ], [ [ 12 ], [ 1/2 ] ] ], 
  [ [ [ 22 ], [ -1/2 ] ], [ [ 19 ], [ 1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], 
      [ [ 24 ], [ -1/2 ] ], [ [ 17 ], [ 1/2 ] ], [ [ 20 ], [ 1/2 ] ], [ [ 21 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 6 ], [ 1/2 ] ], [ [ 3 ], [ -1/2 ] ], 
      [ [ 2 ], [ 1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 1 ], [ -1/2 ] ], 
      [ [ 4 ], [ -1/2 ] ], [ [ 5 ], [ -1/2 ] ], [ [ 13 ], [ 1/2 ] ], [ [  ], [  ] ], 
      [ [ 13 ], [ 1/2 ] ] ], [ [ [ 21 ], [ 1/2 ] ], [ [ 23 ], [ -1/2 ] ], [ [ 20 ], [ 1/2 ] ], 
      [ [ 19 ], [ -1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 18 ], [ 1/2 ] ], 
      [ [ 22 ], [ -1/2 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 5 ], [ -1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 4 ], [ -1/2 ] ], [ [ 3 ], [ 1/2 ] ], 
      [ [ 1 ], [ 1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 2 ], [ -1/2 ] ], [ [ 6 ], [ -1/2 ] ], 
      [ [ 14 ], [ 1/2 ] ], [ [  ], [  ] ], [ [ 14 ], [ 1/2 ] ] ], 
  [ [ [ 19 ], [ 1/2 ] ], [ [ 22 ], [ 1/2 ] ], [ [ 17 ], [ -1/2 ] ], [ [ 21 ], [ 1/2 ] ], 
      [ [ 20 ], [ -1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 24 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [ 3 ], [ -1/2 ] ], 
      [ [ 6 ], [ -1/2 ] ], [ [ 1 ], [ 1/2 ] ], [ [ 5 ], [ -1/2 ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 2 ], [ 1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 7 ], [ -1/2 ] ], [ [ 15 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [ 15 ], [ 1/2 ] ] ], 
  [ [ [ 17 ], [ -1/2 ] ], [ [ 18 ], [ -1/2 ] ], [ [ 19 ], [ -1/2 ] ], [ [ 20 ], [ -1/2 ] ], 
      [ [ 21 ], [ -1/2 ] ], [ [ 22 ], [ -1/2 ] ], [ [ 23 ], [ -1/2 ] ], [ [ 24 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 27 ], [ 1, 1 ] ], [ [ 1 ], [ -1/2 ] ], 
      [ [ 2 ], [ -1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 4 ], [ -1/2 ] ], [ [ 5 ], [ -1/2 ] ], 
      [ [ 6 ], [ -1/2 ] ], [ [ 7 ], [ -1/2 ] ], [ [ 8 ], [ 1/2 ] ], [ [ 16 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [ 16 ], [ 1/2 ] ] ], 
  [ [ [ 16 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 10 ], [ 1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 13 ], [ 1/2 ] ], [ [ 11 ], [ 1/2 ] ], [ [ 9 ], [ -1/2 ] ], 
      [ [ 8 ], [ -1/2 ] ], [ [ 4 ], [ 1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 2 ], [ -1/2 ] ], 
      [ [ 6 ], [ 1/2 ] ], [ [ 5 ], [ -1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 1 ], [ -1/2 ] ], 
      [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 17 ], [ 1/2 ] ], [ [ 17 ], [ 1/2 ] ], 
      [ [  ], [  ] ] ], [ [ [ 12 ], [ 1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 9 ], [ -1/2 ] ], [ [ 11 ], [ 1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 14 ], [ 1/2 ] ], 
      [ [ 10 ], [ -1/2 ] ], [ [ 4 ], [ -1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 5 ], [ 1/2 ] ], 
      [ [ 1 ], [ 1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 6 ], [ -1/2 ] ], 
      [ [ 2 ], [ -1/2 ] ], [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 18 ], [ 1/2 ] ], [ [ 18 ], [ 1/2 ] ], [ [  ], [  ] ] ], 
  [ [ [ 15 ], [ 1/2 ] ], [ [ 13 ], [ 1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 14 ], [ -1/2 ] ], 
      [ [ 10 ], [ -1/2 ] ], [ [ 12 ], [ 1/2 ] ], [ [ 9 ], [ -1/2 ] ], [ [ 11 ], [ -1/2 ] ], 
      [ [ 7 ], [ -1/2 ] ], [ [ 5 ], [ -1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 6 ], [ 1/2 ] ], 
      [ [ 2 ], [ 1/2 ] ], [ [ 4 ], [ -1/2 ] ], [ [ 1 ], [ 1/2 ] ], [ [ 3 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 19 ], [ 1/2 ] ], [ [ 19 ], [ 1/2 ] ], 
      [ [  ], [  ] ] ], [ [ [ 10 ], [ -1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [ 14 ], [ 1/2 ] ], 
      [ [ 16 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 13 ], [ 1/2 ] ], 
      [ [ 12 ], [ -1/2 ] ], [ [ 2 ], [ 1/2 ] ], [ [ 1 ], [ -1/2 ] ], [ [ 6 ], [ -1/2 ] ], 
      [ [ 8 ], [ -1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 3 ], [ 1/2 ] ], [ [ 5 ], [ -1/2 ] ], 
      [ [ 4 ], [ -1/2 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 20 ], [ 1/2 ] ], [ [ 20 ], [ 1/2 ] ], [ [  ], [  ] ] ], 
  [ [ [ 14 ], [ 1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 15 ], [ 1/2 ] ], 
      [ [ 16 ], [ -1/2 ] ], [ [ 9 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 6 ], [ -1/2 ] ], [ [ 3 ], [ 1/2 ] ], [ [ 2 ], [ -1/2 ] ], [ [ 7 ], [ -1/2 ] ], 
      [ [ 8 ], [ -1/2 ] ], [ [ 1 ], [ 1/2 ] ], [ [ 4 ], [ 1/2 ] ], [ [ 5 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 21 ], [ 1/2 ] ], [ [ 21 ], [ 1/2 ] ], 
      [ [  ], [  ] ] ], [ [ [ 13 ], [ -1/2 ] ], [ [ 15 ], [ 1/2 ] ], [ [ 12 ], [ -1/2 ] ], 
      [ [ 11 ], [ 1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 10 ], [ -1/2 ] ], 
      [ [ 14 ], [ -1/2 ] ], [ [ 5 ], [ 1/2 ] ], [ [ 7 ], [ -1/2 ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 3 ], [ -1/2 ] ], [ [ 1 ], [ -1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 2 ], [ 1/2 ] ], 
      [ [ 6 ], [ -1/2 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 22 ], [ 1/2 ] ], [ [ 22 ], [ 1/2 ] ], [ [  ], [  ] ] ], 
  [ [ [ 11 ], [ -1/2 ] ], [ [ 14 ], [ -1/2 ] ], [ [ 9 ], [ 1/2 ] ], [ [ 13 ], [ -1/2 ] ], 
      [ [ 12 ], [ 1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 16 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], 
      [ [ 3 ], [ 1/2 ] ], [ [ 6 ], [ 1/2 ] ], [ [ 1 ], [ -1/2 ] ], [ [ 5 ], [ 1/2 ] ], 
      [ [ 4 ], [ -1/2 ] ], [ [ 2 ], [ -1/2 ] ], [ [ 8 ], [ -1/2 ] ], [ [ 7 ], [ -1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], [ [  ], [  ] ], [ [ 23 ], [ 1/2 ] ], 
      [ [ 23 ], [ 1/2 ] ], [ [  ], [  ] ] ], 
  [ [ [ 9 ], [ -1/2 ] ], [ [ 10 ], [ -1/2 ] ], [ [ 11 ], [ -1/2 ] ], [ [ 12 ], [ -1/2 ] ], 
      [ [ 13 ], [ -1/2 ] ], [ [ 14 ], [ -1/2 ] ], [ [ 15 ], [ -1/2 ] ], [ [ 16 ], [ 1/2 ] ], 
      [ [ 1 ], [ -1/2 ] ], [ [ 2 ], [ -1/2 ] ], [ [ 3 ], [ -1/2 ] ], [ [ 4 ], [ -1/2 ] ], 
      [ [ 5 ], [ -1/2 ] ], [ [ 6 ], [ -1/2 ] ], [ [ 7 ], [ -1/2 ] ], [ [ 8 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [ 25, 26 ], [ 1, 1 ] ], [ [ 24 ], [ 1/2 ] ], 
      [ [ 24 ], [ 1/2 ] ], [ [  ], [  ] ] ], 
  [ [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 9 ], [ 1/2 ] ], [ [ 10 ], [ 1/2 ] ], 
      [ [ 11 ], [ 1/2 ] ], [ [ 12 ], [ 1/2 ] ], [ [ 13 ], [ 1/2 ] ], [ [ 14 ], [ 1/2 ] ], 
      [ [ 15 ], [ 1/2 ] ], [ [ 16 ], [ 1/2 ] ], [ [ 17 ], [ 1/2 ] ], [ [ 18 ], [ 1/2 ] ], 
      [ [ 19 ], [ 1/2 ] ], [ [ 20 ], [ 1/2 ] ], [ [ 21 ], [ 1/2 ] ], [ [ 22 ], [ 1/2 ] ], 
      [ [ 23 ], [ 1/2 ] ], [ [ 24 ], [ 1/2 ] ], [ [ 25 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ] 
     ], [ [ [ 1 ], [ 1/2 ] ], [ [ 2 ], [ 1/2 ] ], [ [ 3 ], [ 1/2 ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 5 ], [ 1/2 ] ], [ [ 6 ], [ 1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 8 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 17 ], [ 1/2 ] ], [ [ 18 ], [ 1/2 ] ], 
      [ [ 19 ], [ 1/2 ] ], [ [ 20 ], [ 1/2 ] ], [ [ 21 ], [ 1/2 ] ], [ [ 22 ], [ 1/2 ] ], 
      [ [ 23 ], [ 1/2 ] ], [ [ 24 ], [ 1/2 ] ], [ [  ], [  ] ], [ [ 26 ], [ 1 ] ], [ [  ], [  ] ] 
     ], [ [ [ 1 ], [ 1/2 ] ], [ [ 2 ], [ 1/2 ] ], [ [ 3 ], [ 1/2 ] ], [ [ 4 ], [ 1/2 ] ], 
      [ [ 5 ], [ 1/2 ] ], [ [ 6 ], [ 1/2 ] ], [ [ 7 ], [ 1/2 ] ], [ [ 8 ], [ 1/2 ] ], 
      [ [ 9 ], [ 1/2 ] ], [ [ 10 ], [ 1/2 ] ], [ [ 11 ], [ 1/2 ] ], [ [ 12 ], [ 1/2 ] ], 
      [ [ 13 ], [ 1/2 ] ], [ [ 14 ], [ 1/2 ] ], [ [ 15 ], [ 1/2 ] ], [ [ 16 ], [ 1/2 ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], 
      [ [ 27 ], [ 1 ] ] ], 1, Zero(F) ]; 
    alg := AlgebraByStructureConstantsArg([F, T, "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
    "j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", 
    "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", 
    "ei", "ej", "ek"], IsSCAlgebraObj and IsJordanAlgebraObj );
    SetFilterObj( alg, IsJordanAlgebra );
    SetJordanRank( alg, 3 );
    SetJordanDegree( alg, 8 );
    # SetJordanMatrixBasis( result.algebra, result.basis );
    SetJordanOffDiagonalBasis( alg, Basis(Oct) );
    SetJordanHomotopeVector( alg, One(alg) );
    SetJordanBasisTraces( alg, List(Basis(alg), 
        j -> (Rank(j)/Dimension(FamilyObj(j)!.fullSCAlgebra))*Trace(AdjointMatrix(Basis(FamilyObj(j)!.fullSCAlgebra), j)))
        );
    jordan_basis := HermitianJordanAlgebraBasis(3, CanonicalBasis(Oct) );
    e := jordan_basis{[1..3]};
    i := jordan_basis{[20..27]};
    j := ComplexConjugate(jordan_basis{[12..19]});
    k := jordan_basis{[4..11]};
    jordan_basis := Concatenation([i,j,k,e]);
    SetJordanMatrixBasis( alg, jordan_basis );
    return alg;
end );

InstallValue( Alb, AlbertAlgebra(Rationals) );

# V := SimpleEuclideanJordanAlgebra(3,8,Basis(Oct));
# B := Basis(V, Concatenation([Basis(V){[20..27]}, -Basis(V){[12..18]}, Basis(V){[19]}, B
# asis(V){[4..11]}, Basis(V){[1..3]}]));
# StructureConstantsTable(B){[1..27]} = StructureConstantsTable(Basis(Alb)){[1..27]};

# T-Design Tools

InstallGlobalFunction( JacobiPolynomial, function(k, a, b)
    local temp, a1, a2, a3, a4, n, x;
    # if not (a > -1 and b > -1) then return fail; fi;
    n := k -1;
    a1 := 2*(n+1)*(n+a+b+1)*(2*n+a+b );
    a2 := (2*n+a+b+1)*(a^2 - b^2 );
    a3 := (2*n + a + b)*(2*n + a + b + 1)*(2*n + a + b + 2 );
    a4 := 2*(n+a)*(n+b)*(2*n+a+b+2 );
    x := [0,1];
    if k = 0 then 
        return [1];
    elif k = -1 then 
        return [0];
    elif k = 1 then 
        return (1/2)*[a-b,(a+b+2)];
    fi;
    return ProductCoeffs([a2/a1, a3/a1], JacobiPolynomial(k-1,a,b))-(a4/a1)*JacobiPolynomial(k-2,a,b );
end );

InstallGlobalFunction( Q_k_epsilon, function(k, epsilon, rank, degree, x) 
    local temp, p, poch, N, m;
    if k = 0 and epsilon = 0 then return 1+0*x; fi;
    m := degree/2;
    N := rank*m;
    p := ValuePol(JacobiPolynomial(k, N-m-1, m-1+epsilon), 2*x - 1 );
    poch := function(a,n)
        if n = 1 then return a;
        elif n =  0 then return 1;
        elif n < 0 then return 0;
        fi;
        return Product(List([1..n], b -> a + b -1) );
    end;
    temp := poch(N, k + epsilon - 1)*poch(N-m, k)*(2*k + N + epsilon - 1 );
    temp := temp/(Factorial(k)*poch(m,k+epsilon) );
    return temp*p/Value(p, [x], [1] );
end );

InstallGlobalFunction( R_k_epsilon, function(k, epsilon, rank, degree, x)
    local temp, n;
    temp := 0*x;
    for n in [0..k] do 
        temp := temp + Q_k_epsilon(n, epsilon, rank, degree, x );
    od;
    return temp;
end );

InstallGlobalFunction( DesignByJordanParameters, function(rank, degree)
    local obj, F, x;
    # Check that the inputs match a Jordan primitive idempotent space    
    if not IsInt(rank) or rank < 2 then return fail; fi;
    if rank > 2 and not degree in [1,2,4,8] then return fail; fi;
    if degree = 8 and rank > 3 then return fail; fi; 
    # Create the object
    obj := Objectify(NewType(NewFamily( "Design"), IsDesign and IsComponentObjectRep, rec()), rec() );
    SetFilterObj(obj, IsAttributeStoringRep );
    # Assign rank and degree attributes.
    SetDesignJordanRank(obj, rank );
    SetDesignJordanDegree(obj, degree );
    # Assign Spherical or Projective filters.
    if rank = 2 then 
        SetFilterObj(obj, IsSphericalDesign );
    fi;
    if degree in [1,2,4,8] then 
        SetFilterObj(obj, IsProjectiveDesign );
    fi;
    return obj;
end );

InstallMethod( PrintObj,
    "for a design",
    [ IsDesign ],
    function(x)
    local text;
    Print( "<design with rank ", DesignJordanRank(x), " and degree ", DesignJordanDegree(x), ">" );
   end );

InstallMethod(DesignQPolynomial,
    "Generic method for designs",
    [ IsDesign ],
    function(D)
        local x, temp;
        x := Indeterminate(Rationals, "x" );
        temp := function(k)
            if IsInt(k) and k > -1 then 
                return CoefficientsOfUnivariatePolynomial((Q_k_epsilon(k, 0, DesignJordanRank(D), DesignJordanDegree(D), x)) );
            fi;
            return fail;
        end;
        return temp;
    end );

InstallMethod( DesignConnectionCoefficients,
    "Generic method for designs",
    [ IsDesign ],
    function(D)
        local temp;
        temp := function(s)
            local x, Q, V, basis, mat, k, i;
            x := Indeterminate(Rationals, "x" );
            Q := List([0..s], i -> Q_k_epsilon(i,0, DesignJordanRank(D), DesignJordanDegree(D), x) );
            V := VectorSpace(Rationals, Q );
            basis := Basis(V, Q );
            mat := [];
            for k in [0..s] do 
                Append(mat, [
                    Coefficients(basis, x^k)
                ] );
            od;
            return mat;
        end;
        return temp; 
    end );

InstallMethod( DesignAddAngleSet, 
    "for designs",
    [ IsDesign, IsList ],
    function(D, A)
        SetDesignAngleSet(D, Set(A) );
        SetFilterObj(D, IsDesignWithAngleSet );
        # Assign Positive Indicator Coefficients filter if applicable.
        if [true] = Set(DesignNormalizedIndicatorCoefficients(D), x -> x > 0) then 
            SetFilterObj(D, IsDesignWithPositiveIndicatorCoefficients );
        fi;
        return D; 
    end );

InstallGlobalFunction( DesignByAngleSet,  function(rank, degree, A)
    local obj, F, x;
    obj := DesignByJordanParameters(rank, degree );
    # Assign angle set.
    SetDesignAngleSet(obj, Set(A) );
    SetFilterObj(obj, IsDesignWithAngleSet ); 
    # Assign Positive Indicator Coefficients filter if applicable.
    if [true] = Set(DesignNormalizedIndicatorCoefficients(obj), x -> x > 0) then 
        SetFilterObj(obj, IsDesignWithPositiveIndicatorCoefficients );
    fi;
    return obj;
end );

InstallMethod( PrintObj,
    "for a design with angle set",
    [ IsDesignWithAngleSet ],
    function(x)
    local text;
    Print( "<design with rank ", DesignJordanRank(x), ", degree ", DesignJordanDegree(x), ", and angle set ", DesignAngleSet(x), ">" );
   end );

InstallMethod( DesignNormalizedAnnihilatorPolynomial,
    "generic method for designs",
    [ IsDesignWithAngleSet ],
    function(D)
    local x, F, A;
    A := DesignAngleSet(D );
    x := Indeterminate(Field(A), "x" );
    F := Product(List(A, a -> (x - a)/(1-a)) );
    return CoefficientsOfUnivariatePolynomial(F );
    end );

InstallMethod( DesignNormalizedIndicatorCoefficients,
    "generic method for designs",
    [ IsDesignWithAngleSet ],
    function(D)
    local x, r, d, Q, F, V, basis;
    r := DesignJordanRank(D );
    d := DesignJordanDegree(D );
    x := Indeterminate(Field(DesignAngleSet(D)), "x" );
    Q := k -> Q_k_epsilon(k, 0, r, d, x );
    F := ValuePol(DesignNormalizedAnnihilatorPolynomial(D), x );
    V := VectorSpace(Field(DesignAngleSet(D)), List([0..Degree(F)], k -> Q(k)) );
    basis := Basis(V, List([0..Degree(F)], k -> Q(k)) );  
    return Coefficients(basis, F );
    end );

InstallMethod( DesignSpecialBound,
    "generic method for designs",
    [ IsDesignWithAngleSet and IsDesignWithPositiveIndicatorCoefficients ],
    function(D)
    if Filtered(DesignNormalizedIndicatorCoefficients(D), x -> x < 0) <> [] then 
        return fail; 
    fi;
    return 1/DesignNormalizedIndicatorCoefficients(D)[1];
    end );

InstallMethod( DesignAddCardinality, 
    "for designs with angle sets",
    [ IsDesignWithAngleSet, IsInt ],
    function(D, v)
        local obj;
        SetDesignCardinality(D, v );
        SetFilterObj(D, IsDesignWithCardinality );
        if DesignSpecialBound(D) = DesignCardinality(D) then 
            SetFilterObj(D, IsSpecialBoundDesign );
            DesignStrength(D );
        fi;
        return D; 
    end );

InstallMethod( PrintObj,
    "for a design with cardinality",
    [ IsDesignWithCardinality ],
    function(x)
        Print( "<design with rank ", DesignJordanRank(x), ", degree ", DesignJordanDegree(x), ", cardinality ", DesignCardinality(x), ", and angle set ", DesignAngleSet(x), ">" );
   end );

InstallMethod( DesignStrength, 
    "method for designs with positive indicator coefficients",
    [ IsDesignWithPositiveIndicatorCoefficients and IsDesignWithCardinality and IsSpecialBoundDesign ],
    function(D)
        local s, i, t, e;
        s := Size(DesignAngleSet(D) );
        for i in [0..s] do 
            if DesignIndicatorCoefficients(D)[i+1] = 1 then 
                t := s + i;
            fi;
        od;
        SetFilterObj(D, IsDesignWithStrength );
        if 0 in DesignAngleSet(D) then 
            e := 1;
        else 
            e := 0;
        fi;
        # Check for tightness, etc
        if t = 2*s - e then
            SetFilterObj(D, IsTightDesign );
            return t;
        fi;
        if t >= 2*s - 2 then
            SetFilterObj(D, IsAssociationSchemeDesign );
            return t;
        fi;
        return t;
    end );

InstallMethod( DesignAnnihilatorPolynomial, 
    "generic method for designs",
    [ IsDesignWithAngleSet and IsDesignWithCardinality ],
    function(D)
        return DesignCardinality(D)*DesignNormalizedAnnihilatorPolynomial(D );
    end );

InstallMethod( DesignIndicatorCoefficients, 
    "generic method for designs",
    [ IsDesignWithAngleSet and IsDesignWithCardinality ],
    function(D)
        return DesignCardinality(D)*DesignNormalizedIndicatorCoefficients(D );
    end );

InstallMethod( PrintObj,
    "for a design with angle set and strength",
    [IsDesignWithAngleSet and IsDesignWithStrength],
    function(x)
    local text;
    Print( "<", DesignStrength(x), "-design with rank ", DesignJordanRank(x), ", degree ", DesignJordanDegree(x), ", cardinality ", DesignCardinality(x), ", and angle set ", DesignAngleSet(x), ">" );
   end );

InstallMethod( PrintObj,
    "for a tight design",
    [IsTightDesign],
    function(x)
        Print( "<Tight ", DesignStrength(x), "-design with rank ", DesignJordanRank(x), ", degree ", DesignJordanDegree(x), ", cardinality ", DesignCardinality(x), ", and angle set ", DesignAngleSet(x), ">" );
   end );

InstallMethod( DesignSubdegrees, 
    "method for a regular scheme design",
    [ IsRegularSchemeDesign and IsDesignWithCardinality and IsDesignWithAngleSet ],
    function(D)
        local rank, degree, v, A, f, s, mat, vec, i;
        v := DesignCardinality(D );
        A := DesignAngleSet(D );
        rank := DesignJordanRank(D );
        degree := DesignJordanDegree(D );
        s := Size(A );
        f := DesignConnectionCoefficients(D)(s );
        # f := ConnectionCoefficients(rank, degree, s );
        mat := [];
        vec := [];
        for i in [0..s-1] do;
            Append(mat, [
                List(Set(A), a -> a^i)
            ] );
            Append(vec, [v*f[i+1][1] - 1] );
        od;
        return SolutionMat(TransposedMat(mat), vec );
    end );

InstallMethod( DesignIntersectionNumbers, 
    "method for an association scheme design",
    [ IsAssociationSchemeDesign ],
    function(D)
        local rank, degree, v, A, F, s, delta, mat, vec, i, j, temp, p, result, gamma, Convolution;
        v := DesignCardinality(D );
        A := DesignAngleSet(D );
        rank := DesignJordanRank(D );
        degree := DesignJordanDegree(D );
        s := Size(A );
        Convolution := function(rank, degree, i, j, x)
            local f, temp, q, y;
            y := Indeterminate(Rationals, "x" );
            q := DesignQPolynomial(D );
            q := k -> ValuePol(DesignQPolynomial(D)(k), y );
            f := DesignConnectionCoefficients(D)(Maximum(i,j) );
            # f := ConnectionCoefficients(rank, degree, Maximum(i,j) );
            temp := List([0..Minimum(i,j)], k -> f[i+1][k+1]*f[j+1][k+1]*q(k) );
            temp := Sum(temp );
            temp := ValuePol(CoefficientsOfUnivariatePolynomial(temp), x );
            return temp;
        end;
        result := [];
        for gamma in A do 
            F := {i,j} -> Convolution(rank, degree, i, j, gamma );
            mat := [];
            vec := [];
            if gamma = 1 then delta := 1; else delta := 0; fi;
            for i in [0..s-1] do 
                for j in [0..s-1] do 
                    Append(mat, [
                        List(Tuples(A,2), a -> a[1]^i *a[2]^j)
                    ] );
                    Append(vec, [v*F(i,j) - gamma^i - gamma^j + delta] );
                    # Append(vec, [v*F(i,j)] );
                od;
            od;
            temp := SolutionMat(TransposedMat(mat), vec );
            p := [];
            for i in [0..s-1] do  
                Append(p, [temp{[1+i*s..s+i*s]}] );
            od;
            Append(result, [p] );
        od;
        for i in [1..s] do 
            result[i] := result[i]+IdentityMat(s+1)*0;
            result[i][i][s+1] := 1;
            result[i][s+1][i] := 1;
        od;
        # degree := DiagonalMat(Subdegrees(rank, degree, v, A)) + IdentityMat(s+1)*0;
        degree := DiagonalMat(DesignSubdegrees(D)) + IdentityMat(s+1)*0;
        degree[s+1][s+1] := 1;
        Append(result, [degree] );
        return result;
    end );

InstallMethod( DesignReducedAdjacencyMatrices, 
    "method for an association scheme design",
    [ IsAssociationSchemeDesign ],
    function(D)
        local s;
        s := Size(DesignAngleSet(D) );
        return List([1..s+1], i -> List([1..s+1], j -> List([1..s+1], k -> DesignIntersectionNumbers(D)[k][i][j])) );
    end );


InstallMethod( DesignBoseMesnerAlgebra, 
    "method for an association scheme design",
    [ IsAssociationSchemeDesign ],
    function(D)
        local p, T, basis, space, coeffs, i, j, k;
        p := DesignReducedAdjacencyMatrices(D );
        space := VectorSpace(Rationals, p );
        basis := Basis(space, p );
        T := EmptySCTable(Length(p), 0, "symmetric" );
        for i in [1..Length(p)] do 
            for j in [1..Length(p)] do 
                coeffs := Coefficients(basis, p[i]*p[j] );
                SetEntrySCTable(T, i, j, Flat(
                    List([1..Length(coeffs)], n -> 
                        [coeffs[n], n]
                    )
                    ) );
            od;
        od;
        return AlgebraByStructureConstants(Rationals, T, "A" );
    end );

InstallMethod( DesignBoseMesnerIdempotentBasis,
    "method for a tight t-design",
    [ IsAssociationSchemeDesign ],
    function(D)
        local i, A, s, v, temp, idempotents, epsilon, final;
        A := DesignAngleSet(D );
        v := DesignCardinality(D );
        s := Size(A );
        if 0 in A then 
            epsilon := 1;
        else 
            epsilon := 0;
        fi;
        idempotents := [];
        # The 0th idempotent
        final := Sum(Basis(DesignBoseMesnerAlgebra(D)))/v;
        # The 1 to s-1 idempotents.
        for i in [1..s-1] do 
            temp := Sum(List([1..s], k -> ValuePol(DesignQPolynomial(D)(i), A[k])*Basis(DesignBoseMesnerAlgebra(D))[k]) ); 
            temp := temp + ValuePol(DesignQPolynomial(D)(i), 1)*Basis(DesignBoseMesnerAlgebra(D))[s + 1];
            temp := temp/v;
            Append(idempotents, [temp] );
        od;
        # The s idempotent.
        temp := One(idempotents) - Sum(idempotents) - final;
        Append(idempotents, [temp] );
        # The 0th idempotent.
        Append(idempotents, [final] );
        idempotents := Basis(DesignBoseMesnerAlgebra(D), idempotents );
        return idempotents;
    end );


InstallMethod( DesignFirstEigenmatrix, 
    "method for a design with association scheme",
    [ IsAssociationSchemeDesign ],
    function(D)
        return List(CanonicalBasis(DesignBoseMesnerAlgebra(D)), x -> Coefficients(DesignBoseMesnerIdempotentBasis(D), x) );
    end );

InstallMethod( DesignSecondEigenmatrix, 
    "method for a design with association scheme",
    [ IsAssociationSchemeDesign ],
    function(D)
        return List(DesignBoseMesnerIdempotentBasis(D), x -> Coefficients(CanonicalBasis(DesignBoseMesnerAlgebra(D)), x))*DesignCardinality(D );
    end );

InstallMethod( DesignMultiplicities, 
    "method for a design with association scheme",
    [ IsAssociationSchemeDesign ],
    function(D)
        local Q, zero; 
        Q := DesignSecondEigenmatrix(D );
        zero := Size(DesignAngleSet(D)) + 1;
        return List([1..zero], i -> Q[i][zero] );
    end );

InstallMethod( DesignValencies, 
    "method for a design with association scheme",
    [ IsAssociationSchemeDesign ],
    function(D)
        local Pmat, zero; 
        Pmat := DesignFirstEigenmatrix(D );
        zero := Size(DesignAngleSet(D)) + 1;
        return List([1..zero], i -> Pmat[i][zero] );
    end );

InstallMethod( DesignKreinNumbers, 
    # Definition in bannai_algebraic_2021 Theorem 2.23, page 61
    "method for a design with association scheme",
    [ IsAssociationSchemeDesign ],
    function(D)
        local s, temp, mat, i, j, k, l, test;
        s := Size(DesignAngleSet(D) );
        temp := [];
        for l in [1..s+1] do
            mat := List([1..s+1], i ->
                List([1..s+1], j -> 
                    (DesignMultiplicities(D)[i]*DesignMultiplicities(D)[j]/DesignCardinality(D))
                    *
                    Sum(
                        List([1..s+1], v ->
                            (1/DesignValencies(D)[v]^2)*
                            DesignFirstEigenmatrix(D)[v][i]*DesignFirstEigenmatrix(D)[v][j]*ComplexConjugate(DesignFirstEigenmatrix(D)[v][l])
                        )
                    )
                )
            );
            Append(temp, [mat] );
        od;
        # Test the result using theorem 2.22(7) on page 59.
        for i in [1..s+1] do 
            for j in [1..s+1] do
                for l in [1..s+1] do 
                    test := DesignSecondEigenmatrix(D)[i][l]*DesignSecondEigenmatrix(D)[j][l] = Sum(List([1..s+1], k -> temp[k][i][j]*DesignSecondEigenmatrix(D)[k][l]) );
                    if test = false then return fail; fi;
                od;
            od;
        od;
        # If test passes, then return the result.
        return temp;     
    end );


# Leech Lattice Tools

InstallGlobalFunction( IsLeechLatticeGramMatrix, function(G)
# Using the classification of integral unimodular lattices, the Leech lattice is the rank 24 unimodular lattice with minimal norm 4.
    local shortest;
    # Confirm M is a basis for a 24 dimensional lattice.
    if not 
        (   
            IsMatrix(G) and 
            DimensionsMat(G) = [24,24] and 
            TransposedMat(G) = G
        )
      then 
        return false;
    fi;
    # Confirm integral lattice (i.e. the lattice is a sublattice of the dual lattice):
    if not Set(Flat(G), x -> IsInt(x)) = [true] then 
        return false;
    fi;
    # Confirm unimodular (i.e. the dual lattice is also a sublattice of the lattice );
    if not Determinant(G) = 1 then 
        return false;
    fi;
    # Confirm no vectors shorter than 4.
    shortest := ShortestVectors(G,3).norms;
    if not shortest = [] then 
        return false;
    fi;
    return true;
end );

InstallGlobalFunction( IsGossetLatticeGramMatrix, function(G)
# Using the classification of integral unimodular lattices, the Gosset lattice is the rank 8 unimodular lattice with minimal norm 2.
    local shortest;
    # Confirm M is a basis for a 8 dimensional lattice.
    if not 
        (   
            IsMatrix(G) and 
            DimensionsMat(G) = [8,8] and 
            TransposedMat(G) = G
        )
      then 
        return false;
    fi;
    # Confirm integral lattice (i.e. the lattice is a sublattice of the dual lattice):
    if not Set(Flat(G), x -> IsInt(x)) = [true] then 
        return false;
    fi;
    # Confirm unimodular (i.e. the dual lattice is also a sublattice of the lattice );
    if not Determinant(G) = 1 then 
        return false;
    fi;
    # Confirm no vectors shorter than 2.
    shortest := ShortestVectors(G,1).norms;
    if not shortest = [] then 
        return false;
    fi;
    return true;
end );

InstallGlobalFunction( OctonionLatticeByGenerators, function(gens, g...)
    local   A,      # The underlying octonion ring.
            obj;    # The resulting lattice object
    # Check that the inputs are correct (a list of equal length lists of octonions from the same algebra implementation in GAP).    
    if 
        not IsOctonionCollColl(gens) or # Needs to be a collection of octonion lists. 
        not IsHomogeneousList(Flat([gens, g])) or # Ensure that the octonions belong to the same octonion algebra. 
        Length(Set(gens, x -> Length(x))) > 1 # The lists need to be equal length 
    then 
        Display( "Usage: OctonionLatticeByGenerators( <gens> [, <g>]) where <gens> is a list of equal length octonion lists and optional argument <g> is a suitable octonion gram matrix." );
        return fail; 
    fi;
    # Construct the Z-module.
    obj := FreeLeftModule(Integers, gens );

    # Determine the octonion algebra.
    A := FamilyObj(One(Flat(gens)))!.fullSCAlgebra;
    SetUnderlyingOctonionRing(obj, A );
    # If no octonion gram matrix is supplied, provide the identity matrix.
    if Length(g) = 0 then 
        SetOctonionGramMatrix(obj, IdentityMat(Length(gens[1]))*One(A) );
    else 
        g := g[1];
        if 
            not IsMatrix(g) or DimensionsMat(gens)[2] <> DimensionsMat(g)[1] 
        then
            Display( "Usage: OctonionLatticeByGenerators( <gens> [, <g>]) where <gens> is a list of equal length octonion lists and optional argument <g> is a suitable octonion gram matrix." );
            return fail;
        fi;
        SetOctonionGramMatrix(obj, g );
    fi;
    # Assign appropriate filters.
    SetFilterObj(obj, IsRowModule );
    SetFilterObj(obj, IsOctonionLattice );
    # Convert the generators into coefficient lists according to the canonical octonion ring basis.
    SetGeneratorsAsCoefficients(obj, 
        List(GeneratorsOfLeftOperatorAdditiveGroup(obj), x -> 
            OctonionToRealVector(CanonicalBasis(UnderlyingOctonionRing(obj)),x)
        )
    );
    # Compute the LLLReducedBasisVectors.
    SetLLLReducedBasisCoefficients(obj, LLLReducedBasis(obj, GeneratorsAsCoefficients(obj)).basis );

    return obj;
end );

InstallMethod( ScalarProduct, 
    "for an octonion lattice",
    [ IsOctonionLattice, IsRowVector, IsRowVector ],
    function(L, x, y)
        local a, b;
        if IsOctonionCollection(x) and IsOctonionCollection(y) then 
            return Trace(x*OctonionGramMatrix(L)*ComplexConjugate(y) );
        fi;
        a := RealToOctonionVector(CanonicalBasis(UnderlyingOctonionRing(L)), x );
        b := RealToOctonionVector(CanonicalBasis(UnderlyingOctonionRing(L)), y );
        return Trace(a*OctonionGramMatrix(L)*ComplexConjugate(b) );
    end );

InstallMethod( GramMatrix,
    "For an octonion lattice", 
    [IsOctonionLattice], 
    function(L)
        return List(LLLReducedBasisCoefficients(L), x ->  List(LLLReducedBasisCoefficients(L), y -> ScalarProduct(L,x,y)) );
    end );

InstallMethod( Rank,
    "For an octonion lattice",
    [IsOctonionLattice],
    function(L)
        return Rank(GeneratorsAsCoefficients(L) );
    end );

InstallMethod( Dimension,
    "For an octonion lattice",
    [IsOctonionLattice],
    function(L)
        return Rank(GeneratorsAsCoefficients(L) );
    end );

InstallMethod( CanonicalBasis,
    "for an octonion lattice",
    [ IsOctonionLattice ],
    function( V )
    local B;
    B:= Objectify( NewType( FamilyObj( V ),
                                IsFiniteBasisDefault
                            and IsCanonicalBasis
                            and IsOctonionLatticeBasis
                            and IsAttributeStoringRep ),
                   rec() );
    SetUnderlyingLeftModule( B, V );
    SetUnderlyingOctonionRing( B, UnderlyingOctonionRing( V ) );
    return B;
    end );

InstallMethod( Basis,
    "for an octonion lattice",
    [ IsOctonionLattice ],
    CanonicalBasis );

InstallMethod( BasisVectors,
    "for canonical basis of a full row module",
    [ IsCanonicalBasis and IsOctonionLatticeBasis ],
    function( B )
        return LLLReducedBasisCoefficients(UnderlyingLeftModule( B ) );
    end );

InstallMethod( TotallyIsotropicCode,
    "For an octonion lattice",
    [ IsOctonionLattice ],
    function(L)
        local lll_basis;
        lll_basis := BasisVectors(CanonicalBasis(L) );
        if Set(Flat(lll_basis), x -> IsInt(x)) = [true] and Set(Flat(GramMatrix(L)*Z(2))) = [Z(2)*0] then 
            return VectorSpace(GF(2), lll_basis*Z(2) );
        fi;
        return fail;
    end );

InstallMethod( \in,
    "for and octonion vector and lattice.",
    IsElmsColls,
    [ IsOctonionCollection, IsOctonionLattice ],
    function( x, L )
        local A;
        A := FamilyObj(One(x))!.fullSCAlgebra;
        if A = UnderlyingOctonionRing(L) then 
            x := OctonionToRealVector(CanonicalBasis(A), x );
            return Set(SolutionMat(LLLReducedBasisCoefficients(L), x), y -> IsInt(y)) = [true];
        fi;
        return false;
    end );

InstallMethod( IsSublattice,
    "For octonion lattices",
    [ IsOctonionLattice, IsOctonionLattice ],
    function(L1, L2)
        if UnderlyingOctonionRing(L1) = UnderlyingOctonionRing(L2) then
            return IsSublattice(LLLReducedBasisCoefficients(L1), LLLReducedBasisCoefficients(L2) );
        fi;
        return false;
    end );

InstallMethod( IsSubset,
    "For octonion lattices",
    IsIdenticalObj,
    [ IsOctonionLattice, IsOctonionLattice ],
        IsSublattice
    );

InstallMethod( \=, 
    "For octonion lattices",
    IsIdenticalObj,
    [ IsOctonionLattice, IsOctonionLattice ],
    function(L1, L2)
        if IsSublattice(L1, L2) and IsSublattice(L2, L1) then 
            return true;
        fi;
        return false;
    end );

# Closure Functions

InstallGlobalFunction( closure_step, function(gens, mult_func, opt...)
    local temp, x, y, z, pair, pairchooser;
    temp := Set(gens, x -> x );
    if Length(opt) > 0 and opt[1] = true then 
        pairchooser := UnorderedTuples;
    else 
        pairchooser := Tuples;
    fi;
    for pair in pairchooser(gens,2) do
        x := pair[1];
        y := pair[2];
        z := mult_func(x,y );
        if not (z in temp) then
            AddSet(temp,z );
        fi;
    od;
    return Set(temp );
end );

InstallGlobalFunction( Closure, function(gens, mult_func, opt...)
    local temp, l;
    if not IsHomogeneousList(gens) then return fail; fi;
    temp := Set(gens );
    l := 0;
    while Length(temp) <> l do
        l := Length(temp );
        if Length(opt) > 0 and opt[1] = true then
            temp := closure_step(temp,mult_func, true );
        else 
            temp := closure_step(temp,mult_func );
        fi;
    od;
    return Set(temp );
end );

InstallGlobalFunction( RandomClosure, function(gens, mult_func, opt...)
    local r, temp, N, n, prior, print_yes;
    temp := ShallowCopy(gens );
    if not IsHomogeneousList(gens) then return fail; fi;
    if Length(opt) > 0 then 
        N := opt[1];
        if Length(opt) > 1 then 
            print_yes := true;
        else 
            print_yes := false;
        fi;
    else 
        N := 1;
        print_yes := false;
    fi;
    n := 0;
    repeat 
        prior := ShallowCopy(temp );
        r := Random(temp );
        temp := Union(temp, List(temp, x -> mult_func(r,x)) );
        if print_yes then 
            Display(Length(temp) );
        fi;
        if temp = prior then 
            n := n+1;
        else 
            n := 0;
        fi;
    until n >= N;
    return temp;
end );

InstallGlobalFunction( RandomOrbitOnSets, function(a_set, start, f, opt...)
    local temp, r, same, limit, print, len;
    if not IsHomogeneousList(Flat([a_set, start])) then return fail; fi;
    if Length(opt)>0 then 
        limit := opt[1];
    else 
        limit := 2;
    fi;
    if Length(opt)>1 then 
        print := true;
    else 
        print := false;
    fi;
    temp := ShallowCopy(start );
    r := Random(a_set );
    same := 0;
    repeat 
        len := Length(temp );
        r := Random(a_set );
        temp := Union(temp, List(temp, x -> f(r,x)) );
        if len = Length(temp) then 
            same := same + 1;
        else 
            same := 0;
        fi;
        if print then 
            Display(len );
        fi;
    until same >= limit;
    return temp;
end );