TSALGEBRA_SteinerTripleSystem_data :=
[ rec(
      autGroup := Group( [ (1,2), (2,3) ] ),
      blocks := [ [ 1, 2, 3 ] ],
      isBinary := true,
      isBlockDesign := true,
      isSTS := true,
      isSimple := true,
      paschConfigCount := 0,
      points := [ 1, 2, 3 ],
      v := 3 ), rec(
      autGroup := Group( [ (1,2)(5,6), (2,3)(6,7), (2,4)(3,5), (4,6)(5,7), (4,5)(6,7) ] ),
      blocks := [ [ 1, 2, 3 ], [ 1, 4, 5 ], [ 1, 6, 7 ], [ 2, 4, 6 ], [ 2, 5, 7 ], [ 3, 4, 7 ], 
          [ 3, 5, 6 ] ],
      isBinary := true,
      isBlockDesign := true,
      isSTS := true,
      isSimple := true,
      paschConfigCount := 7,
      points := [ 1, 2, 3, 4, 5, 6, 7 ],
      v := 7 ), rec(
      autGroup := Group( [ (1,2)(5,6)(7,8), (2,3)(6,9)(7,8), (2,4)(3,5)(8,9), (4,6)(5,7)(8,9), 
          (4,5)(6,8)(7,9) ] ),
      blocks := [ [ 1, 2, 3 ], [ 1, 4, 5 ], [ 1, 6, 7 ], [ 1, 8, 9 ], [ 2, 4, 6 ], [ 2, 5, 8 ], 
          [ 2, 7, 9 ], [ 3, 4, 9 ], [ 3, 5, 7 ], [ 3, 6, 8 ], [ 4, 7, 8 ], [ 5, 6, 9 ] ],
      isBinary := true,
      isBlockDesign := true,
      isSTS := true,
      isSimple := true,
      paschConfigCount := 0,
      points := [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ],
      v := 9 ) ]
;
