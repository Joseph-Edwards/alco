LoadPackage("Alco");

A := OctonionArithmetic(Integers);




# V := SimpleEuclideanJordanAlgebra(3, 8, Basis(A));

# # Work on line 429

S := Subalgebra(Oct, Basis(Oct){[1]});

# comp_alg_basis := Basis(Integers);

# C := UnderlyingLeftModule(comp_alg_basis);
#     F := LeftActingDomain(C);

examples := [
    Basis(Integers),
    Basis(Rationals),
    Basis(CF(3)),
    Basis(AsField(NF(7,[1,6]), CF(7))),
    Basis(AsField(Field(Sqrt(-2)), Field(Sqrt(-2),Sqrt(-5)))),
    Basis(QuaternionAlgebra(Rationals)),
    Basis(QuaternionAlgebra(Integers)),
    Basis(OctonionAlgebra(Rationals)),
    Basis(OctonionAlgebra(Integers)),
    Basis(OctonionArithmetic(Integers)),
    Basis(OctonionArithmetic(Field(Sqrt(5)))),
    Basis(S)
];

J := HermitianSimpleJordanAlgebra(3, Rationals);
for comp_alg_basis in examples do 
    Display(comp_alg_basis);
    herm_mats := HermitianJordanAlgebraBasis(3, comp_alg_basis);
    test1 := List(herm_mats, x -> HermitianMatrixToJordanCoefficients(x, comp_alg_basis)) = IdentityMat(Length(herm_mats));
    Display(test1);
    C := UnderlyingLeftModule(comp_alg_basis);
    K := LeftActingDomain(C);
    # Display([C,K, IsSubset(C, K)]);
    J := HermitianSimpleJordanAlgebra(3, comp_alg_basis, Integers);
    if J = fail then 
        Display("Fail");
    fi;
    test2 := List(herm_mats, x -> HermitianMatrixToJordanVector(x, J)) = BasisVectors(Basis(J));
    Display(test2);
od;
