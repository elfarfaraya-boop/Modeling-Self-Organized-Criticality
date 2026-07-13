/* Avalanche statistics pipeline from CODE SANDPILE.sas — the Size/Time
   avalanche tracking, the call tabulate() frequency tables, the log-log
   transforms, the create/append of the distribution datasets, and the
   PROC REG power-law fit are all exactly as written upstream. The only
   change from the random-driven original is the deposition site: instead of
   a random draw per grain, grains are placed by a fixed deterministic sweep
   so the captured listing is reproducible. Grain count is reduced so the run
   completes quickly. */
proc iml;
   /* parameters */
   n = 50;
   Ngrains = 3000;
   Tas = j(n, n, 0);

   /* storage vectors */
   Size = j(1, Ngrains, 0);
   Time = j(1, Ngrains, 0);

   do i = 1 to Ngrains;

      /* deterministic deposition sweep (replaces the random draw) */
      r = mod( (i-1)*7 , n) + 1;
      c = mod( (i-1)*13, n) + 1;
      Tas[r, c] = Tas[r, c] + 1;

      /* stabilization */
      do while (max(Tas) >= 4);
         idx = loc(Tas >= 4);
         if ncol(idx) = 0 then leave;

      /* variable storage */
      Size[i] = Size[i] + ncol(idx);
      Time[i]  = Time[i] + 1;
         rr = ceil(idx / n);
         cc = idx - (rr - 1)*n;

         /* simultaneous removal */
         do k = 1 to ncol(idx);
            Tas[rr[k], cc[k]] = Tas[rr[k], cc[k]] - 4;
         end;

         /* simultaneous redistribution */
         Add = j(n, n, 0);
         do k = 1 to ncol(idx);
            r = rr[k]; c = cc[k];
            if r > 1  then Add[r-1, c] = Add[r-1, c] + 1;
            if r < n  then Add[r+1, c] = Add[r+1, c] + 1;
            if c > 1  then Add[r, c-1] = Add[r, c-1] + 1;
            if c < n  then Add[r, c+1] = Add[r, c+1] + 1;
         end;

         Tas = Tas + Add;
      end;

   end;

   /* Creation of logarithmic variables and frequency tables */
   Size_pos = Size[loc(Size > 0)]; /* to filter out size = 0  */
   call tabulate(valueS, countS, Size_pos);
   frequencyS = countS/sum(countS);

   log_freq_size = log(countS/sum(countS));
   log_size = log(valueS);

   Time_pos = Time[loc(Time > 0)]; /* to filter out initial Time = 0 */
   call tabulate(valueT, countT, Time_pos);
   frequencyT = countT/sum(countT);

   /* export tables as datasets */

   /* for size */
   create log_size_distribution var {"valueS" "frequencyS" "log_size" "log_freq_size"};
   append;
   close log_size_distribution;

   /* for duration */
   create log_time_distribution var {"valueT" "frequencyT"};
   append;
   close log_time_distribution;

   print (ncol(valueS))[label="n_distinct_sizes"]
         (sum(countS))[label="n_avalanches"]
         (max(valueS))[label="max_avalanche_size"];
quit;

/* Regression: estimate the power-law exponent tau from the log-log fit */
proc reg data=log_size_distribution;
   model log_freq_size = log_size;
run;
quit;
