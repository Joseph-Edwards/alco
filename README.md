# TSALGEBRA

## Setup

1. Install GAP, for Windows or Linux, as you prefer. 
I recommend version 4.9.2 since the latest version seems to have a bug in the installation -  https://www.gap-system.org/Releases/4.9.2.html
 
2. Paste my `tsalgebra` folder (which you’ve downloaded from this private Github project) in your GAP installation as `c:/gap-4.9.2/pkg/tsalgebra`.

3. Open a GAP session and use the command `LoadPackage(“tsalgebra”);`  to import all the commands from this package.

## Instructions

The current version of the package instruction manual is linked here: https://github.com/BNasmith/tsalgebra/blob/master/doc/manual.pdf

## Example Session

```
gap> LoadPackage("tsalgebra");
true
gap> S := SteinerTripleSystem(7,1);;
gap> G := S.autGroup;;
gap> oddSubgroups:=Filtered(ConjugacyClassesSubgroups(G),x->Order(Representative(x)) mod 2 = 1);;
gap> H := Representative(oddSubgroups[Length(oddSubgroups)]);;
gap> IsSubgroup(G,H);
true
gap> STSStarterArcs(S,H);
[ [ 1, 2 ] ]
gap> T := TournamentFromStarterArcs(STSStarterArcs(S,H),H);;
gap> OrientSTS(S,T);
gap> AutGroupOrientedSTS(S);;
gap> Display(S);
rec(
  arcs := [ [ 1, 2 ], [ 1, 5 ], [ 1, 7 ], [ 2, 3 ], [ 2, 4 ], [ 2, 7 ], [ 3, 1 ], [ 3, 4 ],
      [ 3, 5 ], [ 4, 1 ], [ 4, 6 ], [ 4, 7 ], [ 5, 2 ], [ 5, 4 ], [ 5, 6 ], [ 6, 1 ], [ 6, 2 ],
      [ 6, 3 ], [ 7, 3 ], [ 7, 5 ], [ 7, 6 ] ],
  autGroup := Group( [ (1,2,3)(4,5,7), (2,5,7)(3,4,6) ] ),
  blocks := [ [ 1, 2, 3 ], [ 1, 4, 5 ], [ 1, 6, 7 ], [ 2, 4, 6 ], [ 2, 5, 7 ], [ 3, 4, 7 ],
      [ 3, 5, 6 ] ],
  isBinary := true,
  isBlockDesign := true,
  isOrientedSTS := true,
  isSTS := true,
  isSimple := true,
  paschConfigCount := 7,
  points := [ 1 .. 7 ],
  v := 7 )
```
      
