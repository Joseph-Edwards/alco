LoadPackage("Design");

# Convert to integer representations

InstallGlobalFunction(TripleSystemIntegerRepresentation, function(sts)

  local temp, a, pt, block, integer_rep, F, names;
  if IsBlockDesign(sts) then
      F := DefaultField(sts.v);
  else
      F := DefaultField(sts.points);
  fi;
  temp:= rec(points:=[]);

  if F = Rationals then
      return sts;
  else

      if IsPrimeField(F) then
        # ELiminato zero
        integer_rep := function(ffe)
          if ffe = Zero(ffe) then return Length(List(F));
          else return IntFFE(ffe);
          fi;
        end;

      else
        integer_rep := function(ffe)
          if ffe = Zero(ffe) then return Length(List(F));
          else return LogFFE(ffe,PrimitiveElement(F))+1;
          fi;
        end;
      fi;

      # Given our translation tool, integer_rep,
      # Initialize the names translation record
      names := rec();
      # Construct the record
      for pt in sts.points do
          names.(integer_rep(pt)) := pt;
      od;
      temp.pointNames := List([1..Length(RecNames(names))],x->names.(x));

      for pt in sts.points do
        Append(temp.points,[integer_rep(pt)]);
      od;

      # Sort points into set.
      temp.points := Set(temp.points);


      # Relabel the blocks:
      temp.blocks := [];
      for block in sts.blocks do
        Append(temp.blocks,[Set(List(block,x->integer_rep(x)))]);
      od;

      # Sort blocks into set.
      temp.blocks := Set(temp.blocks);

      # Convert point list to a compact range representation and make the point list immutable.
      ConvertToRangeRep(temp.points);
      temp.points := Immutable(temp.points);
      temp.blocks := Immutable(temp.blocks);

      # Include STS order
      temp.v := sts.v;

      if "isSTS" in RecNames(sts) then
          temp.isSTS := sts.isSTS;
      else
          IsSTS(temp);
      fi;

      return temp;
  fi;
end);

InstallGlobalFunction(IsSTS, function(sts)
    local points, blocks, block, pair, a, b, lambda;
    # First check whether sts is a record.
    if not IsRecord(sts) then
        return false;
    fi;
    # Next check whether it has a isSTS status recorded.
    # If it has tested true then return that result.
    if "isSTS" in RecNames(sts) and sts.isSTS then
        return true;
    # Otherwise we'll check for points and blocks to test and test if we can.
    elif "points" in RecNames(sts) and "blocks" in RecNames(sts) then
        # Test the points and blocks for being a Steiner Triple System
        points := sts.points;
        blocks := sts.blocks;
        for pair in Combinations(points,2) do
            lambda := 0;
            a := pair[1]; b:= pair[2];
            for block in blocks do
                if a in block and b in block then
                    lambda := lambda + 1;
                    if lambda > 1 then
                        return false;
                        break;
                    fi;
                fi;
            od;
        od;
        # No false returns by now we can return true.
        sts.isSTS := true;
        # Convert point list to a compact range representation and make the point list immutable.
        ConvertToRangeRep(sts.points);
        sts.points := Immutable(sts.points);
        sts.blocks := Immutable(sts.blocks);

        if not ("v" in RecNames(sts)) then
            sts.v := Length(sts.points);
        fi;
        return true;
    # If we couldn't test we'll return false
    else
        return false;
    fi;
end);

# From Colbourn and Rosa pp. 30-32.
InstallGlobalFunction(PeltesohnTripleSystem, function(n)

    local HDP1, HDP2, B, j, a, b, c, i, m, diff, zero_to_n, temp;

    # Check the input.
    if not IsInt(n) or not (n > 1) or not (n mod 6 in [1,3]) then
        Print("The argument must be a positive integer congruent to 1 or 3 mod 6.\n");
        return fail;
    fi;

    # Define the Heffter difference problem functions

    HDP1 := function(m)

        local result, r, s;

        # Partition the set [1..3*m] into difference triples. The solution is given in Colbourn and Rosa Table 2.3:
        result := [];

        if m = 1 then
            Append(result, [[1,2,3]]);
        elif m = 2 then
            Append(result, [[1,3,4],[2,5,6]]);
        elif m = 3 then
            Append(result, [[1,5,6],[2,8,9],[3,4,7]]);
        elif m mod 3 = 0 then
            s := m / 3;

            for r in [0..s-1] do
                Append(result,[[3*r+1,4*s-r+1,4*s+2*r+2]]);
                Append(result,[[3*r+2,8*s-r,8*s+2*r+2]]);
            od;

            for r in [0..s-2] do
                Append(result,[[3*r+3,6*s-2*r-1,6*s+r+2]]);
            od;

            Append(result,[[3*s,3*s+1,6*s+1]]);

        elif m mod 3 = 1 then
            s := (m - 1)/3;

            for r in [0 .. s-1] do
                Append(result,[[3*r+1,8*s-r+3,8*s+2*r+4]]);
                Append(result,[[3*r+2,6*s-2*r+1,6*s+r+3]]);              Append(result,[[3*r+3,4*s-r+1,4*s+2*r+4]]);
            od;

            Append(result,[[3*s+1,4*s+2,7*s+3]]);

        elif m mod 3 = 2 then
            s := (m - 2)/3;

            for r in [0 .. s-1] do
                Append(result,[[3*r+2,6*s-2*r+3,6*s+r+5]]);
                Append(result,[[3*r+3,8*s-r+5,8*s+2*r+8]]);
            od;

            for r in [0 .. s] do
                Append(result,[[3*r+1,4*s-r+3,4*s+2*r+4]]);
            od;

            Append(result,[[3*s+2,7*s+5,8*s+6]]);
        fi;

        return result;
    end;

    HDP2 := function(m)

        local result, r, s;

        # Partition the set [1..3*m+1]/{2*m + 1} into difference triples. The solution is given in Colbourn and Rosa Table 2.4:
        result := [];

        if m = 2 then
            Append(result, [[1,3,4],[2,6,7]]);
        elif m = 4 then
            Append(result, [[1,12,13],[2,5,7],[3,8,11],[4,6,10]]);
        elif m = 7 then
            Append(result, [[1,11,12],[2,17,19],[3,20,22],[4,10,14],
            [5,8,13],[6,18,21],[7,9,16]]);
        elif m = 10 then
            Append(result, [[1,15,16],[2,27,29],[3,25,28],[4,14,18],
            [5,26,31],[6,17,23],[7,13,20],[8,11,19],[9,24,30],[10,12,22]]);

        elif m mod 3 = 0 then
            s := m / 3;

            for r in [0 .. s-1] do
                Append(result,[[3*r+1,8*s-r+1,8*s+2*r+2]]);
                Append(result,[[3*r+2,4*s-r,4*s+2*r+2]]);
                Append(result,[[3*r+3,6*s-2*r-1,6*s+r+2]]); # SIC there is an error in the book here. 6*s+2*r+2
            od;

        elif m mod 3 = 1 and (m - 1) / 3 > 3 then
            s := (m - 1) / 3;

            for r in [0 .. s] do
                Append(result,[[3*r+1,4*s-r+3,4*s+2*r+4]]);
            od;

            for r in [2 .. s-2] do
                Append(result,[[3*r+2,8*s-r+2,8*s+2*r+4]]);
            od;

            for r in [1 .. s-2] do
                Append(result,[[3*r+3,6*s-2*r+1,6*s+r+4]]);
            od;

            Append(result,[[2,8*s+3,8*s+5],[3,8*s+1,8*s+4],[5,8*s+2,8*s+7],
            [3*s-1,3*s+2,6*s+1],[3*s,7*s+3,8*s+6]]);

        elif m mod 3 = 2 and (m - 2)/3 > 0 then
            s := (m - 2) / 3;

            for r in [0 .. s] do
                Append(result,[[3*r+1,4*s-r+3,4*s+2*r+4]]);
                Append(result,[[3*r+2,8*s-r+6,8*s+2*r+8]]);
            od;

            for r in [0 .. s-1] do
                Append(result,[[3*r+3,6*s-2*r+3,6*s+r+6]]);
            od;

        fi;

        return result;
    end;

    # Determine whether HDP1 or HDP2 is needed and construct the blocks:

    B := [];

    if n mod 6 = 1 then
        m := (n - 1) / 6;

        for diff in HDP1(m) do
            a := diff[1]; b:= diff[2]; c:= diff[3];
            for j in [1..n] do
                Append(B,[[j mod n,(a+j) mod n, (a+b+j) mod n]]);
            od;
        od;

    elif n mod 6 = 3 then
        m := (n - 3) / 6;

        for diff in HDP2(m) do
            a := diff[1]; b:= diff[2]; c:= diff[3];
            for j in [1 .. n] do
                Append(B,[[j mod n,(a+j) mod n, (a+b+j) mod n]]);
            od;
        od;

        # Add the short orbit
        for j in [1 .. n] do
            Append(B, [[j mod n, (2*m+1+j) mod n, (4*m+2+j) mod n]]);
        od;

    fi;

    # Convert the zero point to the n point.
    zero_to_n := function(x)
        if x = 0 then
            return n;
        else
            return x;
        fi;
    end;

    B := List(B, x->List(x,y->zero_to_n(y)));

    # Convert to a set of sets.

    B := Set(List(B,x->Set(x)));

    temp := rec(points:=[1..n],blocks:=B, v:=n, isSTS:=true);

    ConvertToRangeRep(temp.points);
    temp.points := Immutable(temp.points);
    temp.blocks := Immutable(temp.blocks);

    return temp;

end);


InstallGlobalFunction(ProjectiveTripleSystem, function(n)

    local temp, new_block,a,b;

    if not IsInt(n) or n < 1 then
        Print("Requires an integer argument greater than 0.\n");
        return fail;
    fi;

    # Initialize record and set points as elements of a finite field
    temp:= rec(points:=List(GF(2^(n+1))));
    # Remove the Zero.
    Remove(temp.points,1);

    # Initialize list of blocks
    temp.blocks := [];
    for a in temp.points do
        for b in temp.points do
            if a <> b then
                new_block:= Set([a,b,-a-b]);
                    if not (new_block in temp.blocks) then
                    Append(temp.blocks,[new_block]);
                fi;
            fi;
        od;
    od;

    temp.v := Length(temp.points);
    # Certify that this is a triple system based on theory.
    temp.isSTS := true;
# return temp;
    return TripleSystemIntegerRepresentation(temp);

end);


InstallGlobalFunction(AffineTripleSystem, function(n)

  local temp, new_block,a,b;

  if not IsInt(n) or n < 1 then
      Print("Requires an integer argument greater than 0.\n");
      return fail;
  fi;

  # Initialize record and set points as elements of a finite field
  temp:= rec(points:=List(GF(3^n)));

  # Initialize list of blocks
  temp.blocks := [];
  for a in temp.points do
    for b in temp.points do
      if a <> b then
        new_block:= Set([a,b,-a-b]);
        if not (new_block in temp.blocks) then
          Append(temp.blocks,[new_block]);
        fi;
      fi;
    od;
  od;

  temp.v := Length(temp.points);

  # Certify that this is a triple system based on theory.
  temp.isSTS := true;

  # return temp;
  return TripleSystemIntegerRepresentation(temp);
end);

# Retrieve a Netto Triple System data from the library

ReadPackage( "TSALGEBRA", "lib/nettodata.g");

# Retrieve or construct a Netto Triple System

InstallGlobalFunction(NettoTripleSystem, function(q)

  local stored_system, temp, new_block, prime_powers, p, alpha, w, starter, squares,s,b;

  if not IsPrimePowerInt(q) or not (q mod 12 = 7) then
      Print("Requires an prime power integer argument congruent to 7 mod 12.\n");
      return fail;
  fi;

  # Check whether the Netto Triple System is already computed and stored.

  stored_system := Filtered(TSALGEBRA_NettoTripleSystem_data,x->x.v = q);

  # If a system is stored, then stored_system will be a length 1 list. Otherwise it will be length zero.

  if Length(stored_system) = 1 then
      return stored_system[1];
  fi;

  # Find base and exponent:
  prime_powers := PrimePowersInt(q);
  p := prime_powers[1];
  alpha := prime_powers[2];
  # Initialize record and set points as elements of a finite field
  temp:= rec(points:=List(GF(q)));

  # Find the starter block of cube roots.
  w := PrimitiveRoot(GF(q))^((q-1)/3);
  starter := [One(GF(q)),w,w^2];
  temp.starter := starter;

  # Compute the squares:
  squares:=Set(List(GF(q),x -> x^2));
  # Remove the zero.
  Remove(squares,1);

  # Compute blocks
  temp.blocks := [];
  for s in squares do
    for b in temp.points do
      new_block := Set(List(starter,y -> y*s+b));
      if not (new_block in temp.blocks) then
        Append(temp.blocks, [new_block]);
      fi;
    od;
  od;

  temp.v := Length(temp.points);

  # Certify that this is a triple system based on theory.
  temp.isSTS := true;

  # return temp;
  return TripleSystemIntegerRepresentation(temp);
end);

# Retrieve a Steiner Triple System from the library

ReadPackage( "TSALGEBRA", "lib/smallstsdata_with_aut_pasch.g");

# Return all Steiner Triple Systems in the library of order v.
InstallGlobalFunction(SteinerTripleSystems, function(v)
    local temp, n;

    if not IsInt(v) or v < 0 then
        Print("Requires a positive integer argument congruent to 1,3 mod 6.\n");
        return fail;
    fi;


    if (v mod 6) in [1,3] and v > 1 then
        temp := Filtered(TSALGEBRA_SteinerTripleSystem_data,x->x.v = v);
        if Length(temp) = 0 then
            Print("This database does not contain a STS(",String(v),").\n");
        else
            # Certify that this is a triple system based on theory.
            for n in [1..Length(temp)] do
                temp[n].isSTS := true;
                ConvertToRangeRep(temp[n].points);
                temp[n].points := Immutable(temp[n].points);
                temp[n].blocks := Immutable(temp[n].blocks);
            od;
        fi;
        return temp;
    else
        Print("Inadmissible order for a Steiner Triple System.\n");
        return fail;
    fi;
end);

# Return the nth Steiner Triple System in the database of order v.
InstallGlobalFunction(SteinerTripleSystem, function(v,n)
    local temp, sts;
    temp := SteinerTripleSystems(v);
    if temp = fail then
        return fail;
    elif n in [1..Length(temp)] then
        sts := temp[n];
        # Certify that this is a triple system:
        if IsSTS(sts) then
            return sts;
        fi;
    else;
        if Length(temp) = 1 then
            Print("There is only 1 Steiner Triple Systems of order ", String(v), " available in this database.\n");
        elif Length(temp) <> 0 then
            Print("There are ", String(Length(temp))," Steiner Triple Systems of order ", String(v), " available in this database.\n");
        fi;
        return fail;
    fi;
end);


ReadPackage( "TSALGEBRA", "lib/hallHdata.g");

InstallGlobalFunction(HallHTripleSystem, function(n)

    local stored_system, temp, points, blocks, pair, CMLProduct, pointNames;

    if not IsInt(n) or n < 2 then
        Print("Requires an integer argument greater than 2.\n");
        return fail;
    fi;

    # Check whether the Hall H Triple System is already computed and stored.

    stored_system := Filtered(TSALGEBRA_HallH_data,x->x.v = 3^(n+1));

    # If a system is stored, then stored_system will be a length 1 list. Otherwise it will be length zero.

    if Length(stored_system) = 1 then
        return stored_system[1];
    fi;


    CMLProduct := function(x,y)

        local temp, n;

        n := Length(x) - 1;

        temp := List([1..n+1],x->Zero(GF(3)));

        temp[1]:= (x[2]-y[2])*(x[3]*y[4]-x[4]*y[3]);

        return - x - y + temp;

    end;

    points := Tuples(AsList(GF(3)),n+1);
    # Display(points);
    blocks := [];

    for pair in Combinations(points,2) do
        # Display(pair);
        Append(blocks,[[pair[1],pair[2],CMLProduct(pair[1],pair[2])]]);
    od;

    # Convert to integers
    blocks := Set(
                List(blocks,
                    x-> Set(
                            List(x,
                                y->Position(points,y)
                                )
                            )
                    )
                );

    temp := rec();

    pointNames := points;
    points := [1..Length(points)];

    temp.blocks := blocks;
    temp.points := points;
    temp.pointNames := pointNames;
    temp.v := Length(points);

    # Certify that this is a triple system based on theory.
    temp.isSTS := true;

    return temp;
end);

# Triple System Direct Product Contruction: Colbourn and Rosa p. 39

InstallGlobalFunction(STSDirectProduct, function(tsVB,tsWD)

    local V,W,B,D, newpoints, newblocks, point, block, pair, w, relabels, temp;

    if not IsSTS(tsVB) or not IsSTS(tsWD) then
        Print("Requires that both arguments are Steiner Triple Systems.\n");
        return fail;
    fi;

    V:= tsVB.points;
    W:= tsWD.points;
    B:= tsVB.blocks;
    D:= tsWD.blocks;

    newblocks := [];

    # For each {x,y,z} in B, and each i in W, let {x_i,y_i,z_i} be a triple.
    for block in B do
        for point in W do
            Append(newblocks,[List(block,x->[x,point])]);
        od;
    od;

    # For each {i,j,k} in D, i < j < k, and all elements x,y in V (not necessarily distinct), let {x_i,y_j, (x * y)_k} be a triple.
    for block in D do
        block := Set(block);
        for pair in IteratorOfTuples(V,2) do
            w := STSQuasigroupProduct(pair[1],pair[2],tsVB);
            Append(newblocks,[
            [[pair[1],block[1]],
            [pair[2],block[2]],
            [w,block[3]]
            ]]);
        od;

    od;

    # Re-label the blocks.

    newpoints := Set(Concatenation(newblocks));

    relabels := function(y)
        local reference;
        reference := newpoints;
        return Position(reference,y);
    end;

    newblocks := Set(List(newblocks, b-> Set(List(b,p->relabels(p)))));
    newpoints:= List(newpoints,x->relabels(x));

    temp := rec(blocks:=newblocks,points:=newpoints,v:=Length(newpoints));

    # Certify that this is a triple system based on theory.
    temp.isSTS := true;

    return temp;

end);


# Triple System Doubling Contruction: Colbourn and Rosa pp. 41-42 (obscure)
# Clearer in brown1985quadratic 787-788.

InstallGlobalFunction(STSDoublingConstruction, function(sts)

    local temp, V, B, v, newpoints, newblocks, block, x, a, b, c;

    if not IsSTS(sts) then
        Print("The argument must be a Steiner triple system.\n");
        return fail;
    fi;

    V:= sts.points;
    B:= sts.blocks;
    v:= sts.v;

    newpoints := [1..2*v+1];
    newblocks := [];
    # Include Old Blocks
    Append(newblocks,B);


    # Include new blocks joining the two copies:
    for x in V do;
        Append(newblocks,[Set([x,v+1,x+v+1])]);
    od;

    # Final block with two elements in the old and one in the new.
    for block in B do
        a:=block[1];
        b:=block[2];
        c:=block[3];
        Append(newblocks,[Set([a,b+v+1,c+v+1])]);
        Append(newblocks,[Set([b,c+v+1,a+v+1])]);
        Append(newblocks,[Set([c,a+v+1,b+v+1])]);
    od;

    temp := rec(blocks:=newblocks,points:=newpoints,v:=Length(newpoints));
    # Certify that this is a triple system based on theory.
    temp.isSTS := true;
    return temp;

end);

InstallGlobalFunction(STSQuasigroupProduct, function(a,b,sts)

    local block, product, temp;

    # Test the inputs
    if not IsSTS(sts) then
        Print("Third argument must be a Steiner triple system.\n");
        return fail;
    elif not IsInt(a) or not IsInt(b) then
        Print("First two arguments must be integers representing points in a Steiner triple system.\n");
        return fail;
    elif not (a in sts.points) or not (b in sts.points) then
        Print("First two arguments must be integers must be in ", String(sts.points), ".\n");
        return fail;
    fi;

    if a = b then return a; fi;

    temp := Filtered(sts.blocks, x -> (a in x) and (b in x));
    if Length(temp) = 1 then
        block := temp[1];
        product := Filtered(block, x-> not (x in [a,b]))[1];
        return product;
    else
        Print("Error: Check Triple System index. Should be 1.");
        return fail;
    fi;

end);

InstallGlobalFunction(STSLoopProduct, function(a,b,sts)

    local block, product, temp, infty;

    # Test the inputs
    if not IsSTS(sts) then
        Print("Third argument must be a Steiner triple system.\n");
        return fail;
    elif not IsInt(a) or not IsInt(b) then
        Print("First two arguments must be integers representing points in a Steiner loop.\n");
        return fail;
    elif not (a in [1..sts.v+1]) or not (b in [1..sts.v+1]) then
        Print("First two arguments must be integers must be between 1 and ", String(sts.v+1), ".\n");
        return fail;
    fi;

    infty := sts.v + 1;

    if a = b then return infty; fi;

    if a = infty then
        return b;
    elif b = infty then
        return a;
    elif a = b then
        return infty;
    else
        return STSQuasigroupProduct(a,b,sts);
    fi;
end);



InstallGlobalFunction(STSNeighbourhood, function(point,sts)

    local block, pairs, n;

    if not IsSTS(sts) then
        Print("Second argument must be a Steiner Triple System.\n");
        return fail;
    elif not IsInt(point) or not (point in [1..sts.v]) then
        Print("First argument must be a point in ", String(sts.points), ".\n");
        return fail;
    fi;

    pairs := [];

    for block in sts.blocks do
        if point in block then
            Append(pairs,[Set(Filtered(block,x-> not (x=point)))]);
        fi;
    od;

    # If an oriented triple system then we want arcs, not sets.
    if "isOrientedSTS" in RecNames(sts) and sts.isOrientedSTS = true then
        for n in [1..Length(pairs)] do
            if not (pairs[n] in sts.arcs) then
                pairs[n] := Reversed(pairs[n]);
            fi;
        od;
    fi;

    return pairs;

end);

InstallGlobalFunction(STSIncidenceGraph, function(sts)
    local names, edges, blocklabel, block, point, graph, Aut;

    if not IsSTS(sts) then
        Print("The argument must be a Steiner triple system.\n");
        return fail;
    fi;

    names := Concatenation(sts.points,sts.blocks);
    edges := [];

    # Assign edges between points and blocks.
    for block in sts.blocks do
        blocklabel := Position(names,block);
        for point in block do
            Append(edges,[[blocklabel,point]]); #,[point,blocklabel]
            Append(edges,[[blocklabel,point],[point,blocklabel]]);
        od;
    od;

    # Construct Graph
    graph := EdgeOrbitsGraph(Group(()),edges,Length(names));

    # Assign Vertex Names
    AssignVertexNames(graph,names);

    # Find Aut Group and reproduce graph
    Aut := AutGroupGraph(graph);
    StructureDescription(Aut);
    graph := NewGroupGraph(Aut,graph);
    graph.v := sts.v;
    if "pointNames" in RecNames(sts) then
        graph.pointNames:= sts.pointNames;
    fi;

    return graph;
end);


InstallGlobalFunction(STSGraph, function(sts)
    local names, edges, blocklabel, block, point, graphmain, Aut;

    if not IsSTS(sts) then
        Print("The argument must be a Steiner triple system.\n");
        return fail;
    fi;

    names := Concatenation(sts.points,sts.blocks);
    edges := [];

    # Assign edges between points and blocks.
    for block in sts.blocks do
        blocklabel := Position(names,block);
        for point in block do
            Append(edges,[[blocklabel,point]]); #,[point,blocklabel]
            Append(edges,[[blocklabel,point],[point,blocklabel]]);
        od;
    od;

    # Construct Graph
    graphmain := rec(graph:= EdgeOrbitsGraph(Group(()),edges,Length(names)),colourClasses := [[1 .. sts.v],[sts.v+1 .. sts.v+Length(sts.blocks)]]);

    # Assign Vertex Names
    # AssignVertexNames(graphmain.graph,names);

    # Find Aut Group and reproduce graph
    Aut := AutGroupGraph(graphmain);
    StructureDescription(Aut);
    graphmain := rec(graph:=NewGroupGraph(Aut,graphmain.graph),colourClasses :=graphmain.colourClasses);
    graphmain.graph.v := sts.v;

    if "pointNames" in RecNames(sts) then
        graphmain.graph.pointNames:= sts.pointNames;
    fi;

    return graphmain;
end);


# Discussed in Kraski 2004.
STSBlockGraph := function(sts)

    local block_names, pair, a, b, c, temp, edges, graph, Aut;

    if not IsSTS(sts) then
        Print("The argument must be a Steiner triple system.\n");
        return fail;
    fi;

    block_names := sts.blocks;
    edges := [];

    for pair in Combinations(sts.blocks,2) do
        a := ShallowCopy(pair[1]); b:= ShallowCopy(pair[2]);
        temp := ShallowCopy(a);
        IntersectSet(a,b);
        if Length(a) > 0 then
            a := ShallowCopy(temp);
            Append(edges,[[Position(sts.blocks,a),Position(sts.blocks,b)],[Position(sts.blocks,b),Position(sts.blocks,a)]]);
        fi;
    od;
    edges := Set(edges);

    graph := EdgeOrbitsGraph(Group(()),edges, Length(block_names));
    Aut := AutGroupGraph(graph);

    graph := NewGroupGraph(Aut,graph);
    graph.names := block_names;
    StructureDescription(graph.autGroup);
    return graph;

end;

InstallGlobalFunction(STSToBlockDesign, function(sts)

    if not IsSTS(sts) then
        Print("The argument must be a Steiner triple system.\n");
        return fail;
    fi;

    sts.isBlockDesign := true;

end);

InstallGlobalFunction(AutGroupSTS, function(sts)

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    if not IsBlockDesign(sts) then
        STSToBlockDesign(sts);
    fi;

    return AutGroupBlockDesign(sts);

end);

# Return a permutation corresponding to the neighbourhood of a point.

InstallGlobalFunction(STSQuasigroupPermutation,  function(x,sts)

    local temp, nbhd, y, z, pair;

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    if not IsInt(x) or not (x in sts.points) then
        Print("The first argument must be an integer corresponding to a point in a Steiner triple system.\n");
        return fail;
    fi;

    temp := List([1..sts.v]);

    nbhd := STSNeighbourhood(x,sts);

    for pair in nbhd do
        y := pair[1]; z:= pair[2];
        temp[y] := z;
        temp[z] := y;
    od;

    return PermList(temp);
end);

InstallGlobalFunction(STSQuasigroupTranslationGroup, function(sts)

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    return Group(List(sts.points,x->STSQuasigroupPermutation(x,sts)));
end);

InstallGlobalFunction(STSLoopPermutation, function(x,sts)

    local temp, nbhd, y, z, pair;

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    if not IsInt(x) or not (x in [1.. sts.v+1]) then
        Print("The first argument must be an integer corresponding to an element in a Steiner loop.\n");
        return fail;
    fi;

    temp := STSQuasigroupPermutation(x,sts);

    temp := temp * (x,sts.v+1);

    return temp;

end);

InstallGlobalFunction(STSLoopTranslationGroup, function(sts)

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    return Group(List(sts.points,x->STSLoopPermutation(x,sts)));
end);


# Compute the STS point orbits

InstallGlobalFunction(STSPointOrbits, function(sts,H...)

    local G, points, temporbit, results, x;

    # First argument is the triple system.

    if not IsSTS(sts) then
        Print("The first argument of STSBlockOrbits must be a Steiner triple system.\n");
        return fail;
    fi;

    # The optional second argument is a subgroup of the autgroup.

    if "autGroup" in RecNames(sts) then
        G := sts.autGroup;
    else
        G := AutGroupSTS(sts);
    fi;

    if Length(H) = 1 then
        H := H[1];
        if not IsSubgroup(G,H) then
            Print("The optional second argument of STSBlockOrbits must be a subgroup of the Steiner triple system's automorphism group ", G,".\n");
            return fail;
        fi;
    elif Length(H) = 0 then
        H := G;
    else
        Print("The optional second argument of STSBlockOrbits must be a subgroup of the Steiner triple system's automorphism group ", G,".\n");
        return fail;
    fi;

    return Orbits(H,sts.points);

end);



# Compute the STS block orbits

InstallGlobalFunction(STSBlockOrbits, function(sts, H...)

    local G, blocks, orbits, triples, temp, starters, results, x, abc, a, b, c;

    # First argument is the triple system.

    if not IsSTS(sts) then
        Print("The first argument of STSBlockOrbits must be a Steiner triple system.\n");
        return fail;
    fi;

    # The optional second argument is a subgroup of the autgroup.

    if "autGroup" in RecNames(sts) then
        G := sts.autGroup;
    else
        G := AutGroupSTS(sts);
    fi;

    if Length(H) = 1 then
        H := H[1];
        if not IsSubgroup(G,H) then
            Print("The optional second argument of STSBlockOrbits must be a subgroup of the Steiner triple system's automorphism group ", G,".\n");
            return fail;
        fi;
    else
        H := G;
    fi;

    blocks := ShallowCopy(sts.blocks);

    orbits := [];
    triples := [];
    starters := [];

    for abc in sts.blocks do
        a := abc[1]; b := abc[2]; c:= abc[3];
        if not ([a,b,c] in triples) then
            temp := Set(Orbit(H,[a,b,c],OnSets));
            Append(orbits,[temp]);
            Append(triples,temp);
            Append(starters,[temp[1]]);
        fi;
    od;

    return orbits;

end);

InstallGlobalFunction(STSPairOrbits, function(sts,H...)

    local G, orbits, starters, pairs, ab, a, b, temp;

    # First argument is the triple system.

    if not IsSTS(sts) then
        Print("The first argument of STSBlockOrbits must be a Steiner triple system.\n");
        return fail;
    fi;

    # The optional second argument is a subgroup of the autgroup.

    if "autGroup" in RecNames(sts) then
        G := sts.autGroup;
    else
        G := AutGroupSTS(sts);
    fi;

    if Length(H) = 1 then
        H := H[1];
        if not IsSubgroup(G,H) then
            Print("The optional second argument of STSBlockOrbits must be a subgroup of the Steiner triple system's automorphism group ", G,".\n");
            return fail;
        fi;
    else
        H := G;
    fi;

    orbits := [];
    pairs := [];
    starters := [];

    for ab in Iterator(Combinations(sts.points,2)) do
        a := ab[1]; b := ab[2];
        if not ([a,b] in pairs or [b,a] in pairs) then
            temp := Set(Orbit(H,[a,b],OnSets));
            if temp = fail then
                Print("Fails to generate a tournament when taking the orbit of ", [a,b], ".\n");
                return fail;
            else
                Append(orbits,[temp]);
                Append(pairs,temp);
                Append(starters,[temp[1]]);
            fi;
        fi;
    od;

    return orbits;

end);

# Given a permutation group H on set of points V.

InstallGlobalFunction(STSArcOrbits, function(sts,H...)

    local G, orbits, starters, pairs, ab, a, b, temp;

    # First argument is the triple system.

    if not IsSTS(sts) then
        Print("The first argument of STSBlockOrbits must be a Steiner triple system.\n");
        return fail;
    fi;

    # The optional second argument is a subgroup of the autgroup.

    if "autGroup" in RecNames(sts) then
        G := sts.autGroup;
    else
        G := AutGroupSTS(sts);
    fi;

    if Length(H) = 1 then
        H := H[1];
        if not IsSubgroup(G,H) then
            Print("The optional second argument of STSBlockOrbits must be a subgroup of the Steiner triple system's automorphism group ", G,".\n");
            return fail;
        fi;
    else
        H := G;
    fi;

    orbits := [];
    pairs := [];
    starters := [];

    for ab in Iterator(Combinations(sts.points,2)) do
        a := ab[1]; b := ab[2];
        if not ([a,b] in pairs or [b,a] in pairs) then
            temp := Set(Orbit(H,[a,b],OnTuples));
            if temp = fail then
                Print("Fails to generate a tournament when taking the orbit of ", [a,b], ".\n");
                return fail;
            else
                Append(orbits,[temp]);
                Append(pairs,temp);
                Append(starters,[temp[1]]);
            fi;
        fi;
    od;

    return orbits;

end);


# Reverse an orbit in a set of arc orbits.

InstallGlobalFunction(ReverseArcOrbit, function(arc,arcorbits)

    local orbit, n, m, pair;

    if not IsList(arc) or not Length(arc) = 2 then
        Print("The first argument of ReverseArcOrbit must be a representative arc in the orbit.\n");
        return fail;
    fi;

    if not IsList(arcorbits) then
        Print("The second argument of ReverseArcOrbit must be a list of arc orbits.\n");
        return fail;
    else
        for orbit in arcorbits do
            if not Set(List(orbit,x->Length(x))) = Set([2]) then
                Print("The second argument of ReverseArcOrbit must be a list of arc orbits.\n");
                return fail;
            fi;
        od;
    fi;

    for n in [1..Length(arcorbits)] do
        if Reversed(arc) in arcorbits[n] then
            arcorbits[n] := List(arcorbits[n],x->Reversed(x));
        fi;
    od;

end);


# Create a tournament using starter arcs and a permutation group

## NEED TO REVIEW

InstallGlobalFunction(TournamentFromStarterArcs, function(starterarcs,H)

    local arcs, temp;

    if not IsList(starterarcs) or not (Set(List(starterarcs,x->Length(x))) = Set([2])) then
        Print("The first argument of TournamentFromStarterArcs must be a list of starter arcs.\n");
        return fail;
    fi;

    if not IsGroup(H) then
        Print("The second argument of TournamentFromStarterArcs must be a permutation group on the starter arc points.\n");
        return fail;
    fi;

    # Starter arcs are put in set order to ensure that they behave well under the orbit funtion.

    arcs := Set(Concatenation(Orbit(H,Set(starterarcs),OnSetsTuples)));

    temp := Tournament(arcs);

    if temp = fail then
        Print("Failed to form a tournament.");
        return fail;
    else
        temp.starterArcs := arcs;
        temp.group := H;
        return temp;
    fi;
end);


# Create a tournament using starter arcs and a permutation group
InstallGlobalFunction(TournamentFromArcOrbits, function(starterarcs)

    local arcs, temp;

    if not IsList(starterarcs) then
        Print("The first argument of TournamentFromArcOrbits must be a list of arc orbits.\n");
        return fail;
    fi;

    # Starter arcs are put in set order to ensure that they behave well under the orbit funtion.

    arcs := Set(Concatenation(starterarcs));

    temp := Tournament(arcs);

    if temp = fail then
        Print("Failed to form a tournament.");
        return fail;
    else
        return temp;
    fi;
end);

# Given a group, find the odd subgroup conjugacy classes.

InstallGlobalFunction(OddSubgroupsConjugacyClasses,function(G)

    local temp;

    temp := ConjugacyClassesSubgroups(G);
    temp := Filtered(temp, x-> Order(Representative(x)) mod 2 = 1);
    temp := Reversed(temp);

    return temp;

end);

# Given a group, find the odd subgroup conjugacy classes.

InstallGlobalFunction(ReverseArcs,function(arcstoreverse,arcs)

    local temp, arc, n;

    if not IsList(arcstoreverse) or not Set(List(arcstoreverse,x->Length(x))) = [2] or not IsList(arcs) or not Set(List(arcs,x->Length(x))) = [2] then
        Print("The arguments in ReverseArcs must be first a list of arcs to reverse, and second a list of arcs needing correction.\n");
    fi;

    temp := arcstoreverse;

    for arc in temp do
        if arc in arcs then
            n := Position(arcs,arc);
            arcs[n] := Reversed(arcs[n]);
        fi;
    od;

end);


InstallGlobalFunction(STSQuasigroupAlgebra, function(F,sts)

    local T, a, b, c;

    if not IsField(F) then
        Print("The first argument must be a field.\n");
        return fail;
    fi;

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The second argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    STSToBlockDesign(sts);

    # v  elements.
    T := EmptySCTable(sts.v,0,"symmetric");

    for c in [1..sts.v] do
        # For the loop diagonal products c*c = c.
        SetEntrySCTable( T, c, c, [1, c] );
    od;

    # Off-diagonal products:
    for a in [1..sts.v] do
        for b in [1..sts.v] do
            if a <> b then
                c := STSQuasigroupProduct(a,b,sts);
                SetEntrySCTable( T, a, b, [1, c] );
            fi;
        od;
    od;

    return AlgebraByStructureConstants(F,T);

end);

InstallGlobalFunction(STSLoopAlgebra, function(F,sts)

    local T, a, b, c, infty;

    if not IsField(F) then
        Print("The first argument must be a field.\n");
        return fail;
    fi;

    if not IsSTS(sts) and not IsBlockDesign(sts) then
        Print("The second argument must be a Steiner triple system or other block design.\n");
        return fail;
    fi;

    STSToBlockDesign(sts);

    infty := sts.v + 1;

    # v+1  elements.
    T := EmptySCTable(infty,0,"symmetric");

    for c in [1..infty] do
        # For the loop diagonal products c*c = infty.
        SetEntrySCTable( T, c, c, [1, infty] );
        # Also, c*infty = c
        SetEntrySCTable( T, c, infty, [1, c] );
        # Also, infty*c = c
        SetEntrySCTable( T, infty, c, [1, c] );
    od;
        # infty * infty = infty
        SetEntrySCTable( T, infty, infty, [1, infty] );

    # Off-diagonal products:
    for a in [1..sts.v] do
        for b in [1..sts.v] do
            if a <> b then
                c := STSQuasigroupProduct(a,b,sts);
                SetEntrySCTable( T, a, b, [1, c] );
            fi;
        od;
    od;

    return AlgebraByStructureConstants(F,T);

end);

# Create a Left Multiplication Algebra
InstallGlobalFunction(STSLeftAdjointAlgebra, function(F,sts_algebra)

    local A, B, b, temp;

    if not IsField(F) then
        Print("The first argument must be a field.\n");
        return fail;
    fi;

    if not IsAlgebra(sts_algebra) then
        Print("The argument must be an algebra.\n");
        return fail;
    fi;

    B := Basis(sts_algebra);
    temp:=Algebra(Rationals,AdjointBasis(B));

    return temp;

end);

# Generate a tournament record from a list of arcs, or a list or arc orbits.

InstallGlobalFunction(Tournament, function(arcs)

    local temp, T, points, arc, pair, x, y, result;


    # Check that the arcs argument is either a list or arcs or a list or arc orbits.
    if not IsList(arcs) then
        Print("The argument must be a list of arcs or a list or arc orbits.\n");
        return fail;
    elif Length(arcs) = 0 then
        return rec( arcs := [ ], isTournament := true, points := [ 1 ], v := 1 );
    elif Set(List(arcs,x->Length(x))) = [2] then
        T := ShallowCopy(arcs);
    elif Set(List(Concatenation(arcs),x->Length(x))) = [2] then
        T := ShallowCopy(Concatenation(arcs));
    else
        Print("The argument must be a list of arcs or a list or arc orbits.\n");
        return fail;
    fi;

    temp := rec();

    points := Set(Flat(T));

    temp.points := points;
    temp.arcs := T;

    result := true;

    for pair in Combinations(points,2) do
        x := pair[1]; y:= pair[2];

        if [x,y] in T then
            if [y,x] in T then
                result := false;

            fi;
        elif not ([y,x] in T) then
            result := false;

        fi;
    od;

    if result then
        temp.isTournament := true;
        temp.v := Length(temp.points);
        return temp;
    else
        Print("Not a tournament.\n");
        return fail;
    fi;



end);

InstallGlobalFunction(IsTSATournament, function(tournament)

    local T;

    if IsRecord(tournament) and "isTournament" in RecNames(tournament) then
        return tournament.isTournament;
    elif IsRecord(tournament) and "arcs" in RecNames(tournament) then
        T := Tournament(tournament.arcs);
        if T = fail then
            return false;
        else
            tournament.isTournament := T.isTournament;
            return T.isTournament;
        fi;
    else
        T := Tournament(tournament);
        if T = fail then
            return false;
        else
            return T.isTournament;
        fi;
    fi;

end);


InstallGlobalFunction(RandomTSATournament, function(v)

    local arcs, pair, x, y, coin;

    arcs := [];

    for pair in Combinations([1..v],2) do
        x := pair[1]; y := pair[2];

        coin := Random(GlobalRandomSource,0,1);

        if coin = 0 then
            Append(arcs,[[x,y]]);
        else
            Append(arcs,[[y,x]]);
        fi;
    od;

    return Tournament(arcs);

end);


InstallGlobalFunction(TournamentToDigraph, function(tournament)

    local graph, aut;

    if not IsTSATournament(tournament) then Print("Requires tournament input"); return fail; fi;

    graph := EdgeOrbitsGraph(Group([()]),tournament.arcs,Maximum(tournament.points));

    aut := AutGroupGraph(graph);
    StructureDescription(aut);

    graph:= NewGroupGraph(aut,graph);

    graph.v := tournament.v;

    return graph;

end);

InstallGlobalFunction(AutGroupTournament, function(tournament)

    local digraph;
    digraph := TournamentToDigraph(tournament);
    tournament.autGroup := digraph.autGroup;

    return digraph.autGroup;

end);

InstallGlobalFunction(OrientSTS, function(sts,tournament)

    sts.isOrientedSTS := true;
    sts.arcs := tournament.arcs;

end);

InstallGlobalFunction(IsOrientedSTS, function(sts)

    if not ("isOrientedSTS" in RecNames(sts)) then
        return false;
    else
        return sts.isOrientedSTS;
    fi;

end);

InstallGlobalFunction(RemoveSTSOrientation, function(sts)

    local temp;

    if IsOrientedSTS(sts) then
        Unbind(sts.arcs);
        Unbind(sts.isOrientedSTS);
    fi;

end);

InstallGlobalFunction(OrientedSTSIncidenceGraph, function(sts)

    local triv, I, edges, temp, names;

    if not IsOrientedSTS(sts) then
        Print("The argument must be an oriented Steiner triple system.\n");
        return fail;
    fi;

    triv := Group([()]);

    I := STSIncidenceGraph(sts);

    edges := [];

    Append(edges,DirectedEdges(I));
    Append(edges,sts.arcs);

    names := Concatenation(sts.points,sts.blocks);

    temp := EdgeOrbitsGraph(triv,edges,Length(I.names));

    AssignVertexNames(temp,names);

    return temp;

end);

InstallGlobalFunction(AutGroupOrientedSTS, function(sts)

    local I, temp, gens, aut;

    if not IsOrientedSTS(sts) then
        Print("The argument must be an oriented Steiner triple system.\n");
        return fail;
    fi;

    I := OrientedSTSIncidenceGraph(sts);

    temp := AutGroupGraph(I);

    if Size(temp) = 1 then
        aut := temp;
    else
        gens := List(GeneratorsOfGroup(temp),g->Permutation(g,sts.points));
        aut := Group(gens);
        StructureDescription(aut);
    fi;

    sts.autGroup := aut;

    return aut;

end);

InstallGlobalFunction(STSLoopTwistedAlgebra, function(F, sts)

    local SCT, v, infty, points, blocks, prod, x, y, z, pair, arcs, zeta, zetastar, oriented_graph;

    if not IsField(F) then
        Print("The first argument must be a field.\n");
        return fail;
    fi;

    if not IsOrientedSTS(sts) then
        Print("The second argument must be an oriented Steiner triple system or other block design.\n");
        return fail;
    fi;

    oriented_graph := OrientedSTSIncidenceGraph(sts);

    # Determine algebra dimension
    v := sts.v;
    infty := v + 1;

    # Initialize Structure Constant Table
    SCT := EmptySCTable(infty, 0);

    # List of points:
    points := Filtered(oriented_graph.names,x-> not IsList(x));

    # List of blocks:
    blocks := Filtered(oriented_graph.names,x-> IsList(x));

    # Quasigroup product on blocks:
    prod := function(a,b,blocks)
        local c, block, temp;
        for block in blocks do
            if (a in block) and (b in block) then
                block := Set(block);
                RemoveSet(block,a);
                RemoveSet(block,b);
                c := block[1];
                return c;
            fi;
        od;
    end;

    # Set the factor system for infty and symmetric products:
    SetEntrySCTable( SCT, infty, infty, [1, infty] );

    for x in points do
        SetEntrySCTable( SCT, x, x, [-1, infty] );
        SetEntrySCTable( SCT, x, infty, [1, x] );
        SetEntrySCTable( SCT, infty, x, [1, x] );
    od;

    # Set the factor system for arcs:
    arcs := DirectedEdges(InducedSubgraph(oriented_graph,points));

    zeta := 1;
    zetastar := -1;

    for pair in arcs do
        x := pair[1]; y:= pair[2];
        z := prod(x, y, blocks);
        SetEntrySCTable( SCT, x, y, [zeta, z]);
        SetEntrySCTable( SCT, y, x, [zetastar, z]);
    od;



    return AlgebraByStructureConstants(F,SCT);

end);

InstallGlobalFunction(STSQuasigroupTwistedAlgebra, function(F, sts)

    local SCT, v, infty, points, blocks, prod, x, y, z, pair, arcs, zeta, zetastar, oriented_graph;

    if not IsField(F) then
        Print("The first argument must be a field.\n");
        return fail;
    fi;

    if not IsOrientedSTS(sts) then
        Print("The second argument must be an oriented Steiner triple system or other block design.\n");
        return fail;
    fi;

    oriented_graph := OrientedSTSIncidenceGraph(sts);

    # Determine algebra dimension
    v := sts.v;

    # Initialize Structure Constant Table
    SCT := EmptySCTable(v, 0);

    # List of points:
    points := Filtered(oriented_graph.names,x-> not IsList(x));

    # List of blocks:
    blocks := Filtered(oriented_graph.names,x-> IsList(x));

    # Quasigroup product on blocks:
    prod := function(a,b,blocks)
        local c, block, temp;
        for block in blocks do
            if (a in block) and (b in block) then
                block := Set(block);
                RemoveSet(block,a);
                RemoveSet(block,b);
                c := block[1];
                return c;
            fi;
        od;
    end;

    # Set the factor system for infty and symmetric products:
    for x in points do
        SetEntrySCTable( SCT, x, x, [-1, x] );
    od;

    # Set the factor system for arcs:
    arcs := DirectedEdges(InducedSubgraph(oriented_graph,points));

    zeta := 1;
    zetastar := -1;

    for pair in arcs do
        x := pair[1]; y:= pair[2];
        z := prod(x, y, blocks);
        SetEntrySCTable( SCT, x, y, [zeta, z]);
        SetEntrySCTable( SCT, y, x, [zetastar, z]);
    od;



    return AlgebraByStructureConstants(F,SCT);

end);

# Generate a Cyclotomic Steiner Loop Algebra from an Oriented STS Incidence Graph
InstallGlobalFunction(STSLoopCyclotomicAlgebra, function(sts, n)

    local infty, SCT, zeta, zetastar, pair, x, y, z;


    if not IsOrientedSTS(sts) then
        Print("The first argument must be an oriented Steiner triple system.\n");
        return fail;
    fi;

    if not IsInt(n) or n < 1 then
        Print("The second argument must be a positive integer.\n");
        return fail;
    fi;
    # Determine algebra dimension
    infty := sts.v + 1;

    # Initialize Structure Constant Table
    SCT := EmptySCTable(infty, 0);

    # Set the factor system for infty and symmetric products:
    SetEntrySCTable( SCT, infty, infty, [1, infty] );

    for x in sts.points do
        SetEntrySCTable( SCT, x, x, [1, infty] );
        SetEntrySCTable( SCT, x, infty, [1, x] );
        SetEntrySCTable( SCT, infty, x, [1, x] );
    od;

    zeta := E(n);
    zetastar := ComplexConjugate(zeta);

    for pair in sts.arcs do
        x := pair[1]; y:= pair[2];
        z := STSLoopProduct(x, y, sts);
        SetEntrySCTable( SCT, x, y, [zeta, z]);
        SetEntrySCTable( SCT, y, x, [zetastar, z]);
    od;

    return AlgebraByStructureConstants(Cyclotomics,SCT);

end);

# Generate a Cyclotomic Steiner Loop Algebra from an Oriented STS Incidence Graph
InstallGlobalFunction(STSQuasigroupCyclotomicAlgebra, function(sts, n)

    local infty, SCT, zeta, zetastar, pair, x, y, z;


    if not IsOrientedSTS(sts) then
        Print("The first argument must be an oriented Steiner triple system.\n");
        return fail;
    fi;

    if not IsInt(n) or n < 1 then
        Print("The second argument must be a positive integer.\n");
        return fail;
    fi;

    # Initialize Structure Constant Table
    SCT := EmptySCTable(sts.v, 0);


    for x in sts.points do
        SetEntrySCTable( SCT, x, x, [1, x] );
    od;

    zeta := E(n);
    zetastar := ComplexConjugate(zeta);

    for pair in sts.arcs do
        x := pair[1]; y:= pair[2];
        z := STSQuasigroupProduct(x, y, sts);
        SetEntrySCTable( SCT, x, y, [zeta, z]);
        SetEntrySCTable( SCT, y, x, [zetastar, z]);
    od;

    return AlgebraByStructureConstants(Cyclotomics,SCT);

end);

InstallGlobalFunction(STSPaschConfigCount, function(sts)

    local pasch, four_blocks, temp, block;

    if not IsSTS(sts) then
        Print("The first argument must be a Steiner triple system.\n");
        return fail;
    fi;

    if "paschConfigCount" in RecNames(sts) then
        return sts.paschConfigCount;
    else

        pasch := 0;

        for four_blocks in Combinations(sts.blocks,4) do
            temp := four_blocks[1];
            for block in four_blocks do
                temp:= UnionSet(temp,block);
            od;
            if Length(temp) = 6 then
                pasch := pasch + 1;
            fi;
        od;

        sts.paschConfigCount := pasch;

        return pasch;

    fi;

end);


# Classify the block-subtournaments as either cyclic or transitive

InstallGlobalFunction(ClassifyBlockSubtournaments, function(sts)

    local block, i, a, b, c, C3, T3, orientation;

    if not IsOrientedSTS(sts) then
        Print("The argument must be an oriented Steiner triple system.\n");
        return fail;
    fi;

    # Will create an block orientation list. Each entry gives us 1 for cyclic standard, -1 for cyclic reversed, 0 for transitive.

    orientation := List(sts.blocks,x->3);

    for block in sts.blocks do
        i := Position(sts.blocks,block);
        a := block[1]; b:= block[2]; c:= block[3];
        if [a,b] in sts.arcs then
            if [b,c] in sts.arcs then
                if [c,a] in sts.arcs then
                    orientation[i] := 1;
                else
                    orientation[i] := 0;
                fi;
            else
                orientation[i] := 0;
            fi;
        else
            if [a,c] in sts.arcs then
                if [c,b] in sts.arcs then
                    orientation[i] := -1;
                else
                    orientation[i] := 0;
                fi;
            else
                orientation[i] := 0;
            fi;
        fi;
    od;

    sts.blockOrientation := orientation;

    return orientation;

end);

# Apply an orientation to the blocks of a triple system. 1 leaves them in set order. -1 reverses them.

InstallGlobalFunction(OrientSTSBlocks, function(sts, list)

    local i, a, b, c;

    if not IsSTS(sts) then
        Print("The first argument must be a Steiner triple system.\n");
        return fail;
    fi;

    if not IsList(list) or Length(list) <> Length(sts.blocks) or not (Set(list) = Set([-1,1]))  then
        Print("The second argument must be list with 1 or -1 entries of length", Length(sts.blocks),".\n");
        return fail;
    fi;

    sts.arcs := [];

    for i in [1..Length(sts.blocks)] do
        a := sts.blocks[i][1]; b:= sts.blocks[i][2]; c:= sts.blocks[i][3];
        if list[i] = 1 then
            Append(sts.arcs, [[a,b],[b,c],[c,a]]);
        elif list[i] = -1 then
            Append(sts.arcs, [[a,c],[c,b],[b,a]]);
        fi;
    od;

    sts.blockOrientation := list;
    sts.isOrientedSTS := true;

end);

# Apply a random cyclic orientation to the blocks of an STS.

InstallGlobalFunction(RandomlyOrientSTSBlocks, function(sts)

    local orientation;

    if not IsSTS(sts) then
        Print("The first argument must be a Steiner triple system.\n");
        return fail;
    fi;

    orientation := List(sts.blocks,x->Random(GlobalRandomSource,[-1,1]));

    OrientSTSBlocks(sts,orientation);

end);

# A prepared Oriented Fano Plane.
InstallValue(OrientedFanoPlane,
rec(
  arcs := [ [ 1, 2 ], [ 1, 3 ], [ 1, 5 ], [ 2, 3 ], [ 2, 4 ], [ 2, 6 ],
      [ 3, 4 ], [ 3, 5 ], [ 3, 7 ], [ 4, 1 ], [ 4, 5 ], [ 4, 6 ],
      [ 5, 2 ], [ 5, 6 ], [ 5, 7 ], [ 6, 1 ], [ 6, 3 ], [ 6, 7 ],
      [ 7, 1 ], [ 7, 2 ], [ 7, 4 ] ], autGroup := Group([ (1,2,3,4,5,6,7), (2,3,5)(4,7,6) ]),
  blocks := [ [ 1, 2, 4 ], [ 1, 3, 7 ], [ 1, 5, 6 ], [ 2, 3, 5 ],
      [ 2, 6, 7 ], [ 3, 4, 6 ], [ 4, 5, 7 ] ], isBinary := true,
  isBlockDesign := true, isOrientedSTS := true, isSTS := true,
  isSimple := true,
  pointNames := [ Z(7)^0, Z(7)^2, Z(7), Z(7)^4, Z(7)^5, Z(7)^3, 0*Z(7) ],
  points := [ 1 .. 7 ], v := 7 )
);

InstallValue(OrientedTriple,
rec(
  arcs := [ [ 1, 2 ], [ 2, 3 ], [ 3, 1 ] ],
  autGroup := Group( [ (1,2), (2,3) ] ),
  blocks := [ [ 1, 2, 3 ] ],
  isBinary := true,
  isBlockDesign := true,
  isOrientedSTS := true,
  isSTS := true,
  isSimple := true,
  paschConfigCount := 0,
  points := [ 1 .. 3 ],
  v := 3 )
);

InstallGlobalFunction(OctonionAlgebra,function(F)

    local temp, SCT;
    # return AlgebraByStructureConstants(F,SCT,
    # "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
    # "j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", 
    # "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", 
    # "e1", "e2", "e3");


    temp := STSLoopTwistedAlgebra(F,OrientedFanoPlane);

    SCT := StructureConstantsTable(Basis(temp));

    return AlgebraByStructureConstants( F, SCT,  
    "i1", "i2", "i3", "i4", "i5", "i6", "i7", "e");

end);


InstallGlobalFunction(LeftMultiplicationMatrix, function(B,v)

    local
        bb,     # Basis vectors
        n,      # Dimension of basis
        F,      # Underlying field
        i,      # Index
        temp;

        if not (IsBasis(B) and IsVector(v)) then
            Print("#W The argument must be an algebra basis and a vector.");
            return fail;
        fi;

        bb := BasisVectors(B);
        n := Length(bb);
        temp := [];

        for i in [1..n] do
            Append(temp,[Coefficients(B,v*bb[i])]);
        od;

        return temp;

end);

InstallGlobalFunction(RightMultiplicationMatrix, function(B,v)

    local
        bb,     # Basis vectors
        n,      # Dimension of basis
        F,      # Underlying field
        i,      # Index
        temp;

        if not (IsBasis(B) and IsVector(v)) then
            Print("#W The argument must be an algebra basis and a vector.");
            return fail;
        fi;

        bb := BasisVectors(B);
        n := Length(bb);
        temp := [];

        for i in [1..n] do
            Append(temp,[Coefficients(B,bb[i]*v)]);
        od;

        return temp;

end);

InstallGlobalFunction(LeftMultiplicationAlgebra, function(B)
    local
        leftgens,   # Matrix representations of the left basis.
        bb,         # Basis vectors of B.
        n,          # Dimension of B.
        F,          # Field corresponding to the basis.
        i,          # Index.
        temp;

    if not IsBasis(B) then
        Print("#W The argument must be an algebra basis.");
        return fail;
    fi;

    bb := BasisVectors(B);
    n := Length(bb);
    F:= LeftActingDomain( UnderlyingLeftModule( B ) );
    leftgens := [];

    # Append the matrix element of each basis vector.
    for i in [1..n] do
        Append(leftgens,[LeftMultiplicationMatrix(B,bb[i])]);
    od;

    return Algebra(F,leftgens);

end);

InstallGlobalFunction(RightMultiplicationAlgebra, function(B)
    local
        rightgens,  # Matrix representations of the right basis.
        bb,         # Basis vectors of B.
        n,          # Dimension of B.
        F,          # Field corresponding to the basis.
        i,          # Index.
        temp;

    if not IsBasis(B) then
        Print("#W The argument must be an algebra basis.");
        return fail;
    fi;

    bb := BasisVectors(B);
    n := Length(bb);
    F:= LeftActingDomain( UnderlyingLeftModule( B ) );
    rightgens := [];

    # Append the matrix element of each basis vector.
    for i in [1..n] do
        Append(rightgens,[RightMultiplicationMatrix(B,bb[i])]);
    od;

    return Algebra(F,rightgens);

end);

InstallGlobalFunction(FullMultiplicationAlgebra, function(B)
    local
        gens,  # Matrix representations of the right basis.
        bb,         # Basis vectors of B.
        n,          # Dimension of B.
        F,          # Field corresponding to the basis.
        i,          # Index.
        temp;

    if not IsBasis(B) then
        Print("#W The argument must be an algebra basis.");
        return fail;
    fi;

    bb := BasisVectors(B);
    n := Length(bb);
    F:= LeftActingDomain( UnderlyingLeftModule( B ) );
    gens := [];

    # Append the matrix element of each basis vector.
    for i in [1..n] do
        Append(gens,[LeftMultiplicationMatrix(B,bb[i])]);
        Append(gens,[RightMultiplicationMatrix(B,bb[i])]);
    od;

    return Algebra(F,Set(gens));

end);

InstallGlobalFunction(NettoQuasigroupProduct, function(u,v)

    local p, alpha, q, epsilon, not_square;

    # Check that u,v are finite field elements of equal characteristic.
    if not (IsFFE(u) and IsFFE(v))
        or Characteristic(u) <> Characteristic(v)
        or Characteristic(u) mod 12 <> 7
        or DegreeFFE(u) mod 2 <> 1
        or DegreeFFE(v) mod 2 <> 1
        then
        Print("Both arguments must be finite field elements of equal characteristic with prime power order congruent to 7 mod 12.\n");
        return fail;
    fi;

    if u = v then
        return u;
    fi;

    # Find the prime power of the system.
    p := Characteristic(u);
    alpha := Maximum([DegreeFFE(u),DegreeFFE(v)]);
    q := p^alpha;

    # Find the sixth roots of unity

    epsilon := [Z(p)^((p-1)/6),Z(p)^((p-1)*5/6)];

    if Length(Filtered(epsilon,x->x^2-x+1 = Zero(Z(q)))) <> 2 then
        Print("Error finding the sixth roots of unity.\n");
    fi;

    # Determine whether u-v is a square.

    not_square := LogFFE(u-v,Z(q)) mod 2;

    if not_square = 0 then
        return u*epsilon[1] + v*epsilon[2];
    elif not_square = 1 then
        return v*epsilon[1] + u*epsilon[2];
    fi;

end);


InstallGlobalFunction(CayleyDicksonTournamentDoubling, function(T)

    local p, infty, new_arcs, a;

    if not IsTSATournament(T) then
        Print("The argument must be a tournament.\n");
        return fail;
    fi;

    infty := T.v + 1;

    new_arcs := [];

    for p in T.points do
        Append(new_arcs,[[p,p+infty]]);
        Append(new_arcs,[[p+infty,infty]]);
        Append(new_arcs,[[infty,p]]);
    od;

    for a in T.arcs do
        Append(new_arcs,[a]);
        Append(new_arcs,[[a[2]+infty,a[1]+infty]]);
        Append(new_arcs,[[a[2],a[1]+infty]]);
        Append(new_arcs,[[a[2]+infty,a[1]]]);
    od;

    return Tournament(Set(new_arcs));

end);


InstallGlobalFunction(CayleyDicksonOrientedSTSDoubling, function(sts)

    local temp;

    if not IsOrientedSTS(sts) then
        Print("The argument must be an oriented STS.\n");
        return fail;
    fi;

    temp := STSDoublingConstruction(sts);

    OrientSTS(temp,CayleyDicksonTournamentDoubling(sts));

    return temp;

end);

InstallGlobalFunction(TournamentInNeighbourhoods, function(T)

    local temp, x;

    if not IsTSATournament(T) then
        Print("The argument must be a tournament.\n");
        return fail;
    fi;

    temp := rec();

    for x in T.points do

        temp.(x) := Set(List(Filtered(T.arcs,arc->arc[2] = x),y -> y[1]));

    od;

    return temp;

end);

InstallGlobalFunction(TournamentOutNeighbourhoods, function(T)

    local temp, x;

    if not IsTSATournament(T) then
        Print("The argument must be a tournament.\n");
        return fail;
    fi;

    temp := rec();

    for x in T.points do

        temp.(x) := Set(List(Filtered(T.arcs,arc->arc[1] = x),y -> y[2]));

    od;

    return temp;

end);


InstallGlobalFunction(TransitiveTournaments, function(T,n)

    local temp, trans, m, t, out, nextpts, next;

    if not IsTSATournament(T) then
        Print("The first argument must be a tournament.\n");
        return fail;
    fi;

    if not IsInt(n) and n > 0 then
        Print("The second argument must be a positive integer.\n");
        return fail;
    fi;

    # Create a record of transitive tournaments of a particular length
    trans := rec();
    trans.(1) := Set(T.points);
    trans.(2) := Set(T.arcs);

    if n < 3 then
        return trans.(n);
    fi;

    out := TournamentOutNeighbourhoods(T);

    m := 3;

    while m < n+1 do

        temp := [];



        for t in trans.(m-1) do

            nextpts := Intersection(List(t,p->out.(p)));

            if Length(nextpts) > 0 then

                for next in nextpts do

                    Append(temp, [Concatenation(t,[next])]);

                od;

            fi;

        od;

        trans.(m) := Set(temp);

        m := m + 1;

    od;

    return trans;

end);


InstallGlobalFunction(TransitiveTournamentCounter, function(T,max)

    local temp, m, listversion;

    if not IsTSATournament(T) then
        Print("The first argument must be a tournament.\n");
        return fail;
    fi;

    if not IsInt(max) and max > 0 then
        Print("The second argument must be a positive integer.\n");
        return fail;
    fi;

    listversion := TransitiveTournaments(T,max);

    temp := rec();

    for m in RecNames(listversion) do

        temp.(m) := Length(listversion.(m));

    od;

    return temp;

end);


# Construct the Albert Algebra from a Stored Structure Constant Table 

ReadPackage( "TSALGEBRA", "lib/albertalgebraSCT.g");

InstallGlobalFunction(AlbertAlgebra, function(F)

    local temp, SCT;
  
    SCT := ALBERTALGEBRAStructureConstants;

    return AlgebraByStructureConstants(F,SCT,
    "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
    "j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", 
    "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", 
    "e1", "e2", "e3");
end);

ReadPackage( "TSALGEBRA", "lib/AlbertAlgebraIntegralRoots.g");

ReadPackage( "TSALGEBRA", "lib/OctonionUnitLoopsAnLRB.g");


ReadPackage( "TSALGEBRA", "lib/staralgebraSCT.g");

InstallGlobalFunction(StarAlgebra, function(F)

    local temp, SCT;
  
    SCT := STARPRODUCTStructureConstants;

    return AlgebraByStructureConstants(F,SCT,
    "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
    "j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", 
    "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8");

    # return AlgebraByStructureConstants(F,StarProductStructureConstants(F), 
    # "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", 
    # "j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", 
    # "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8");
end);

# Commutative star product on 24 dimensions.

InstallGlobalFunction(StarFlatProduct, function(a,b,F)
    # Input two vectors of length 24.
    local temp, x1,x2,x3, y1,y2,y3, z1,z2,z3, Oct;  
    Oct := OctonionAlgebra(F);
    x1 := LinearCombination(Basis(Oct),a{[1..8]});
    x2 := LinearCombination(Basis(Oct),a{[9..16]});
    x3 := LinearCombination(Basis(Oct),a{[17..24]});
    y1 := LinearCombination(Basis(Oct),b{[1..8]});
    y2 := LinearCombination(Basis(Oct),b{[9..16]});
    y3 := LinearCombination(Basis(Oct),b{[17..24]});
    z1 := Coefficients(Basis(Oct),x2*y3 + y2*x3);
    z2 := Coefficients(Basis(Oct),x3*y1 + y3*x1);
    z3 := Coefficients(Basis(Oct),x1*y2 + y1*x2);
    return Concatenation([z1,z2,z3]);
end);

InstallGlobalFunction(StarProductStructureConstants, function(F)
    local T, coeffs, i, j, u, v, temp;
    # Empty Structure Constant Table
    T := EmptySCTable(24, Zero(F), "symmetric");
    for i in [1..24] do 
        for j in [i..24] do 
        # Construct two vectors to multiply
            u := List([1..24],x->Zero(F));
            u[i] := One(F);
            v := List([1..24],x->Zero(F));
            v[j] := One(F);
        # Multiply them.
            temp := StarFlatProduct(u,v,F);
        # Find the coefficents of the product vector.
            coeffs := List([1..24],x->[temp[x],x]);
            coeffs := Flat(coeffs);
        # Add an entry to the SCTable.
            SetEntrySCTable( T, i, j, coeffs );
        od;
    od;
    return T;
end);
