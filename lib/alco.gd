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

# Quaternion Tools

DeclareGlobalVariable( "QuaternionD4Basis" );

# Icosian and Golden Field Tools

DeclareGlobalFunction( "GoldenRationalComponent" );

DeclareGlobalFunction( "GoldenIrrationalComponent" );

DeclareGlobalFunction( "GoldenModSigma" );

DeclareGlobalVariable( "IcosianH4Basis" );

# Octonion Algebra and Arithmetic Tools

DeclareCategory( "IsOctonion", IsScalar );

DeclareCategory( "IsOctonionArithmeticElement", IsOctonion );

DeclareCategoryCollections( "IsOctonion" );

DeclareCategoryCollections( "IsOctonionCollection" );

DeclareAttribute( "Norm", IsOctonionCollection );

BindGlobal( "OctonionAlgebraData", NEW_SORTED_CACHE(true) );

DeclareCategory( "IsOctonionAlgebra", IsOctonionCollection and IsFullSCAlgebra );

    DeclareAttribute( "GramMatrix", IsOctonionAlgebra );

DeclareOperation("IsOctavianInt", [ IsOctonion ] );

DeclareGlobalName( "OctavianIntegers" );

DeclareCategory( "IsOctavianIntegers", IsFLMLOR
  and IsFiniteDimensional );

DeclareRepresentation(
    "IsCanonicalBasisOctavianIntegersRep", IsAttributeStoringRep );

DeclareGlobalFunction( "OctonionAlgebra" );

DeclareGlobalVariable( "OctonionE8Basis" );

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

DeclareGlobalVariable( "Alb" );

DeclareGlobalFunction( "HermitianMatrixToAlbertVector" );

DeclareGlobalFunction( "AlbertVectorToHermitianMatrix" );

DeclareOperation("P", [IsJordanAlgebraObj, IsJordanAlgebraObj]);

DeclareOperation("P", [IsJordanAlgebraObj]);

DeclareOperation("JTS", [IsJordanAlgebraObj, IsJordanAlgebraObj, IsJordanAlgebraObj]);

# T-Design Tools

DeclareGlobalFunction( "JacobiPolynomial" );

DeclareGlobalFunction( "Q_k_epsilon" );

DeclareGlobalFunction( "R_k_epsilon" );

# Designs

DeclareCategory( "IsDesign", IsObject );

DeclareCategory( "IsSphericalDesign", IsDesign );

DeclareCategory( "IsProjectiveDesign", IsDesign );

DeclareGlobalFunction( "DesignByJordanParameters" );

DeclareAttribute( "DesignJordanRank", IsDesign );

DeclareAttribute( "DesignJordanDegree", IsDesign );

DeclareAttribute( "DesignQPolynomial", IsDesign );

DeclareAttribute( "DesignConnectionCoefficients", IsDesign );

DeclareCategory( "IsDesignWithAngleSet", IsDesign );

DeclareAttribute( "DesignAngleSet", IsDesignWithAngleSet );

DeclareOperation( "DesignAddAngleSet", [ IsDesign, IsList ] );

DeclareGlobalFunction( "DesignByAngleSet" );

DeclareAttribute( "DesignNormalizedAnnihilatorPolynomial", IsDesignWithAngleSet );

DeclareAttribute( "DesignNormalizedIndicatorCoefficients", IsDesignWithAngleSet );

DeclareCategory( "IsDesignWithPositiveIndicatorCoefficients", IsDesignWithAngleSet );

DeclareAttribute( "DesignSpecialBound", IsDesignWithAngleSet ); 

DeclareCategory( "IsDesignWithCardinality", IsDesign );

DeclareCategory( "IsRegularSchemeDesign", IsDesignWithCardinality );

DeclareCategory( "IsSpecialBoundDesign", IsRegularSchemeDesign );

DeclareCategory( "IsAssociationSchemeDesign", IsSpecialBoundDesign );

DeclareCategory( "IsTightDesign", IsAssociationSchemeDesign );

DeclareAttribute( "DesignCardinality", IsDesignWithAngleSet );

DeclareOperation( "DesignAddCardinality", [ IsDesignWithAngleSet, IsInt ] );

DeclareCategory( "IsDesignWithStrength", IsDesign );

DeclareAttribute( "DesignAnnihilatorPolynomial", IsDesignWithAngleSet and IsDesignWithCardinality );

DeclareAttribute( "DesignIndicatorCoefficients", IsDesignWithAngleSet and IsDesignWithCardinality );

DeclareAttribute( "DesignStrength", IsDesignWithAngleSet and IsDesignWithCardinality );

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
    
    DeclareAttribute( "OctonionGramMatrix", IsOctonionLattice );

    DeclareAttribute( "GeneratorsAsCoefficients", IsOctonionLattice );

    DeclareAttribute( "LLLReducedBasisCoefficients", IsOctonionLattice );

    DeclareAttribute( "GramMatrix", IsOctonionLattice );

    DeclareAttribute( "TotallyIsotropicCode", IsOctonionLattice );

DeclareGlobalFunction( "OctonionLatticeByGenerators" );

DeclareCategory( "IsOctonionLatticeBasis", IsCanonicalBasis );

    DeclareAttribute( "UnderlyingOctonionRing", IsOctonionLatticeBasis );

DeclareOperation( "IsSublattice", [ IsOctonionLattice, IsOctonionLattice ] );

# Other Tools

DeclareGlobalFunction( "Closure" );

DeclareGlobalFunction( "RandomClosure" );

DeclareGlobalFunction( "RandomOrbitOnSets" );

