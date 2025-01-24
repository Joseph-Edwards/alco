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

# Cyclotomic Tools

DeclareOperation( "IsEisenInt", [ IsCyc ] );

DeclareGlobalName( "EisensteinIntegers");

DeclareCategory( "IsEisensteinIntegers", IsEuclideanRing and IsFLMLOR
  and IsFiniteDimensional );

DeclareRepresentation(
    "IsCanonicalBasisEisensteinIntegersRep", IsAttributeStoringRep,
    [ "conductor", "zumbroichbase" ] );

DeclareOperation( "Coefficients", [ IsCanonicalBasisEisensteinIntegersRep,
      IsCyc ] );

DeclareOperation("IsKleinInt", [ IsCyc ] );

DeclareGlobalName( "KleinianIntegers");

DeclareCategory( "IsKleinianIntegers", IsEuclideanRing and IsFLMLOR
  and IsFiniteDimensional );

DeclareRepresentation(
    "IsCanonicalBasisKleinianIntegersRep", IsAttributeStoringRep );

DeclareOperation( "Coefficients", [ IsCanonicalBasisKleinianIntegersRep,
      IsCyc ] );

# Quaternion Tools

DeclareGlobalVariable( "QuaternionD4Basis" );

DeclareOperation("IsHurwitzInt", [ IsQuaternion ] );

DeclareGlobalName( "HurwitzIntegers" );

DeclareCategory( "IsHurwitzIntegers", IsFLMLOR
  and IsFiniteDimensional and IsQuaternionCollection );

DeclareRepresentation(
    "IsCanonicalBasisHurwitzIntegersRep", IsAttributeStoringRep );

DeclareOperation( "Coefficients", [ IsCanonicalBasisHurwitzIntegersRep,
      IsQuaternion ] );

# Icosian and Golden Field Tools

DeclareGlobalFunction( "GoldenRationalComponent" );

DeclareGlobalFunction( "GoldenIrrationalComponent" );

DeclareGlobalFunction( "GoldenModSigma" );

DeclareGlobalVariable( "IcosianH4Generators" );

DeclareOperation("IsIcosian", [ IsQuaternion ] );

DeclareGlobalName( "IcosianRing" );

DeclareCategory( "IsIcosianRing", IsFLMLOR
  and IsFiniteDimensional );

DeclareRepresentation(
    "IsCanonicalBasisIcosianRingRep", IsAttributeStoringRep );

DeclareOperation( "Coefficients", [ IsCanonicalBasisIcosianRingRep,
      IsQuaternion ] );

# Octonion Algebra and Arithmetic Tools

DeclareCategory( "IsOctonion", IsScalar );

DeclareCategory( "IsOctonionArithmeticElement", IsOctonion );

DeclareCategoryCollections( "IsOctonion" );

DeclareCategoryCollections( "IsOctonionCollection" );

DeclareAttribute( "Norm", IsOctonionCollection );

BindGlobal( "OctonionAlgebraData", NEW_SORTED_CACHE(true) );

DeclareCategory( "IsOctonionAlgebra", IsOctonionCollection and IsFullSCAlgebra );

    DeclareAttribute( "GramMatrix", IsOctonionAlgebra );

DeclareGlobalVariable( "OctonionE8Basis" );

DeclareOperation("IsOctavianInt", [ IsOctonion ] );

DeclareGlobalName( "OctavianIntegers" );

DeclareCategory( "IsOctavianIntegers", IsFLMLOR
  and IsFiniteDimensional );

DeclareRepresentation(
    "IsCanonicalBasisOctavianIntegersRep", IsAttributeStoringRep );

DeclareOperation( "Coefficients", [ IsCanonicalBasisOctavianIntegersRep,
      IsOctonion ] );

DeclareGlobalFunction( "OctonionAlgebra" );

DeclareGlobalFunction( "OctonionToRealVector" );

DeclareGlobalFunction( "RealToOctonionVector" );

DeclareGlobalFunction( "VectorToIdempotentMatrix" );

DeclareGlobalFunction( "WeylReflection" );

# Jordan Algebra Tools

DeclareCategory( "IsJordanAlgebra", IsAlgebra );

    DeclareAttribute( "JordanRank", IsJordanAlgebra );

    DeclareAttribute( "JordanDegree", IsJordanAlgebra );

    DeclareAttribute( "JordanMatrixBasis", IsJordanAlgebra );

    DeclareAttribute( "JordanOffDiagonalBasis", IsJordanAlgebra );

    DeclareAttribute( "JordanHomotopeVector", IsJordanAlgebra );

    DeclareAttribute( "JordanBasisTraces", IsJordanAlgebra );

DeclareCategory( "IsJordanAlgebraObj", IsSCAlgebraObj );

    DeclareGlobalFunction( "GenericMinimalPolynomial" );

    DeclareAttribute( "JordanRank", IsJordanAlgebraObj );

    DeclareAttribute( "JordanDegree", IsJordanAlgebraObj );

    DeclareAttribute( "JordanAdjugate", IsJordanAlgebraObj );

    DeclareAttribute( "IsPositiveDefinite", IsJordanAlgebraObj );

DeclareGlobalFunction( "HermitianJordanAlgebraBasis" );

DeclareGlobalFunction( "HermitianMatrixToJordanCoefficients" );

DeclareGlobalFunction( "HermitianMatrixToJordanVector" );

DeclareGlobalFunction( "JordanSpinFactor" );

DeclareGlobalFunction( "HermitianSimpleJordanAlgebra" );

DeclareAttribute( "JordanAlgebraGramMatrix", IsJordanAlgebra );

DeclareGlobalFunction( "JordanHomotope" );

DeclareGlobalFunction( "SimpleEuclideanJordanAlgebra" );

# Albert Algebra Tools

BindGlobal( "AlbertAlgebraData", NEW_SORTED_CACHE(true) );

DeclareGlobalFunction( "AlbertAlgebra" );

DeclareCategory( "IsAlbertAlgebraObj", IsJordanAlgebraObj );

DeclareGlobalFunction( "HermitianMatrixToAlbertVector" );

DeclareGlobalFunction( "AlbertVectorToHermitianMatrix" );

DeclareOperation("JordanQuadraticOperator", [IsJordanAlgebraObj, IsJordanAlgebraObj]);

DeclareOperation("JordanQuadraticOperator", [IsJordanAlgebraObj]);

DeclareOperation("JordanTripleSystem", [IsJordanAlgebraObj, IsJordanAlgebraObj, IsJordanAlgebraObj]);

# T-Design Tools

DeclareGlobalFunction( "JacobiPolynomial" );

DeclareGlobalFunction( "Q_k_epsilon" );

DeclareGlobalFunction( "R_k_epsilon" );

# Designs

DeclareCategory( "IsJordanDesign", IsObject );

DeclareCategory( "IsSphericalDesign", IsJordanDesign );

DeclareCategory( "IsProjectiveDesign", IsJordanDesign );

DeclareGlobalFunction( "DesignByJordanParameters" );

DeclareAttribute( "DesignJordanRank", IsJordanDesign );

DeclareAttribute( "DesignJordanDegree", IsJordanDesign );

DeclareAttribute( "DesignQPolynomials", IsJordanDesign );

DeclareAttribute( "DesignConnectionCoefficients", IsJordanDesign );

DeclareCategory( "IsJordanDesignWithAngleSet", IsJordanDesign );

DeclareAttribute( "DesignAngleSet", IsJordanDesignWithAngleSet );

DeclareOperation( "DesignAddAngleSet", [ IsJordanDesign, IsList ] );

DeclareGlobalFunction( "DesignByAngleSet" );

DeclareAttribute( "DesignNormalizedAnnihilatorPolynomial", IsJordanDesignWithAngleSet );

DeclareAttribute( "DesignNormalizedIndicatorCoefficients", IsJordanDesignWithAngleSet );

DeclareCategory( "IsJordanDesignWithPositiveIndicatorCoefficients", IsJordanDesignWithAngleSet );

DeclareAttribute( "DesignSpecialBound", IsJordanDesignWithAngleSet ); 

DeclareCategory( "IsJordanDesignWithCardinality", IsJordanDesign );

DeclareCategory( "IsRegularSchemeDesign", IsJordanDesignWithCardinality );

DeclareCategory( "IsSpecialBoundDesign", IsRegularSchemeDesign );

DeclareCategory( "IsAssociationSchemeDesign", IsSpecialBoundDesign );

DeclareCategory( "IsTightDesign", IsAssociationSchemeDesign );

DeclareAttribute( "DesignCardinality", IsJordanDesignWithAngleSet );

DeclareOperation( "DesignAddCardinality", [ IsJordanDesignWithAngleSet, IsInt ] );

DeclareCategory( "IsJordanDesignWithStrength", IsJordanDesign );

DeclareAttribute( "DesignAnnihilatorPolynomial", IsJordanDesignWithAngleSet and IsJordanDesignWithCardinality );

DeclareAttribute( "DesignIndicatorCoefficients", IsJordanDesignWithAngleSet and IsJordanDesignWithCardinality );

DeclareAttribute( "DesignStrength", IsJordanDesignWithAngleSet and IsJordanDesignWithCardinality );

DeclareAttribute( "DesignSubdegrees", IsRegularSchemeDesign );

DeclareAttribute( "DesignIntersectionNumbers", IsAssociationSchemeDesign );

DeclareAttribute( "DesignReducedAdjacencyMatrices", IsAssociationSchemeDesign );

DeclareAttribute( "DesignBoseMesnerAlgebra", IsAssociationSchemeDesign );

DeclareAttribute( "DesignBoseMesnerIdempotentBasis", IsAssociationSchemeDesign );

DeclareAttribute( "DesignFirstEigenmatrix", IsAssociationSchemeDesign );

DeclareAttribute( "DesignSecondEigenmatrix", IsAssociationSchemeDesign );

DeclareAttribute( "DesignMultiplicities", IsAssociationSchemeDesign );

DeclareAttribute( "DesignValencies", IsAssociationSchemeDesign );

DeclareAttribute( "DesignKreinNumbers", IsAssociationSchemeDesign );

# Octonion Lattice Tools

DeclareGlobalFunction( "IsLeechLatticeGramMatrix" );

DeclareGlobalFunction( "IsGossetLatticeGramMatrix" );

DeclareGlobalVariable( "MOGLeechLatticeGeneratorMatrix" );

DeclareGlobalVariable( "MOGLeechLatticeGramMatrix" );

DeclareCategory( "IsOctonionLattice", IsFreeLeftModule );

    DeclareAttribute( "UnderlyingOctonionRing", IsOctonionLattice );

    DeclareAttribute( "UnderlyingOctonionRingBasis", IsOctonionLattice );
    
    DeclareAttribute( "OctonionGramMatrix", IsOctonionLattice );

    DeclareAttribute( "GeneratorsAsCoefficients", IsOctonionLattice );

    DeclareAttribute( "LLLReducedBasisCoefficients", IsOctonionLattice );

    DeclareAttribute( "GramMatrix", IsOctonionLattice );

    DeclareAttribute( "TotallyIsotropicCode", IsOctonionLattice );

DeclareGlobalFunction( "OctonionLatticeByGenerators" );

DeclareCategory( "IsOctonionLatticeBasis", IsCanonicalBasis );

    DeclareAttribute( "UnderlyingOctonionRing", IsOctonionLatticeBasis );

DeclareOperation( "IsSublattice", [ IsOctonionLattice, IsOctonionLattice ] );

DeclareOperation( "Coefficients", [ IsOctonionLatticeBasis,
      IsOctonionCollection ] );

# Other Tools

DeclareGlobalFunction( "Closure" );

DeclareGlobalFunction( "RandomClosure" );

DeclareGlobalFunction( "RandomOrbitOnSets" );

