/* Optimized (center-driven) BTW sandpile — the loc()-vectorized toppling from
   CODE SANDPILE.sas, on a smaller grid/grain count so it completes quickly.
   The loc() cell-selection, the ceil(idx/n) / idx-(rr-1)*n index arithmetic,
   the simultaneous 4-grain removal, and the Add-matrix redistribution are
   exactly as written upstream. */
proc iml;
/* Model parameters */
n = 51;
centre = ceil(n/2);

Tas = j(n, n, 0);
Tas[centre, centre] = 2000;

/* Critical threshold */
do while (max(Tas) >= 4);
   idx = loc(Tas >= 4);
   if ncol(idx)=0 then leave;

   /* Positions (r,c) of these cells */
   rr = ceil(idx / n);
   cc = idx - (rr - 1)*n; *or "cc = mod(idx - 1, n) + 1";

   /* Remove 4 grains simultaneously */
   do k = 1 to ncol(idx);
      Tas[ rr[k], cc[k] ] = Tas[ rr[k], cc[k] ] - 4;
   end;

   /* Accumulate redistributed grains */
   Add = j(n,n,0);

   do k = 1 to ncol(idx);
      r = rr[k];  c = cc[k];
      if r>1  then Add[r-1,c] = Add[r-1,c] + 1;  /* up   */
      if r<n  then Add[r+1,c] = Add[r+1,c] + 1;  /* down */
      if c>1  then Add[r,c-1] = Add[r,c-1] + 1;  /* left */
      if c<n  then Add[r,c+1] = Add[r,c+1] + 1;  /* right */
   end;

   Tas = Tas + Add;
end;

/* Visualisation values: grain counts mod 4 (the fractal colouring) */
Sandpile = mod(Tas, 4);
maxz = max(Tas);
totz = sum(Tas);
print maxz[label="max_grains_any_cell"] totz[label="total_grains"];
print (Sandpile[centre, (centre-3):(centre+3)])[label="center_row_mod4"];
quit;
