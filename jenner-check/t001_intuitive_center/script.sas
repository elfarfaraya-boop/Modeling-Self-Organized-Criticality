/* Intuitive (center-driven) BTW sandpile — the nested-loop toppling from
   CODE SANDPILE.sas, run with a smaller grain count so it completes quickly.
   The threshold rule, the +1 neighbour redistribution, and the do-while
   stabilization are exactly as written upstream. */
proc iml;
/* Create a 50x50 matrix of zeros */
Tas = j(50, 50, 0);
n = nrow(Tas); /* number of rows = 50 */
centre = ceil(n/2);

/* Loop: drop grains one by one */
do Ngrains = 1 to 500;
   Tas[centre, centre] = Tas[centre, centre] + 1;

   /* Topple all cells >= 4 */
   repeat = 1;
   do while (repeat);
      repeat = 0; /* assume no more topplings */

      do r = 1 to 50;
         do c = 1 to 50;

            /* If a cell has 4 grains or more : it topples */
            if Tas[r,c] >= 4 then do;
               Tas[r,c] = Tas[r,c] - 4;

               /* Give +1 to neighbouring cells if they exist */
               if r > 1  then Tas[r-1,c] = Tas[r-1,c] + 1; /* UP */
               if r < 50 then Tas[r+1,c] = Tas[r+1,c] + 1; /* DOWN */
               if c > 1  then Tas[r,c-1] = Tas[r,c-1] + 1; /* LEFT */
               if c < 50 then Tas[r,c+1] = Tas[r,c+1] + 1; /* RIGHT */

               /* Reset repeat=1 since some cells may now topple */
               repeat = 1;
            end;
         end;
      end;
   end;
end;

/* Report the stabilized configuration in numbers (deterministic) */
maxz = max(Tas);
totz = sum(Tas);
print (Ngrains-1)[label="grains_dropped"] maxz[label="max_grains_any_cell"] totz[label="total_grains"];
print (Tas[centre, (centre-3):(centre+3)])[label="center_row_window"];
quit;
