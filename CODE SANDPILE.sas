/*SAS IML implementation of the intuitive version*/
proc iml;
/* Create a 50x50 matrix of zeros */
Tas = j(50, 50, 0);
n = nrow(Tas); /* number of rows = 50 */
centre = ceil(n/2);

/* Loop: drop Ngrains one by one */
do Ngrains = 1 to 4000;
   Tas[centre, centre] = Tas[centre, centre] + 1;

   /* 3. Topple all cells >= 4 */
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

print tas;

call heatmapdisc(Tas)
   xvalues = 1:n
   yvalues = 1:n
   title = "Sandpile 50x50 4000 grains"; 
quit;


/*SAS IML implementation of the optimised version*/

proc iml;
/* Model parameters */
n = 201;
centre = ceil(n/2);

Tas = j(n, n, 0);
Tas[centre, centre] = 400000;             

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

/* Visualisation: Heatmap */
Sandpile = mod(Tas, 4);
call heatmapdisc(Sandpile)
   xvalues=1:n  yvalues=1:n
   title="Abelian sandpile (201x201, 400000 grains)"
   displayoutlines=0;
quit;


/*Random Sandpile Code*/
proc iml;
   /*parameters*/
   n = 50;                            
   Ngrains = 200000;                  
   Tas = j(n, n, 0);                 

   /* storage vectors */
   Size = j(1, Ngrains, 0);           
   Time  = j(1, Ngrains, 0);          

   do i = 1 to Ngrains;

      /* random deposition */
      r = ceil(n * rand("Uniform"));
      c = ceil(n * rand("Uniform"));
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

   /*final visualization*/
   Sandpile = mod(Tas, 4);
   call heatmapdisc(Sandpile)
        xvalues=1:n yvalues=1:n
        title="Abelian Sandpile (random)"
        displayoutlines=0;



/*Creation of logarithmic variables and frequency tables*/
Size_pos = Size[loc(Size > 0)]; /* to filter out size = 0  */
call tabulate(valueS, countS, Size_pos);
frequencyS = countS/sum(countS);

log_freq_size = log(countS/sum(countS));
log_size = log(valueS);

Time_pos = Time[loc(Time > 0)]; /* to filter out initial Time = 0 */
call tabulate(valueT, countT, Time_pos);
frequencyT = countT/sum(countT);

/*export tables as datasets*/

/* for size */
create log_size_distribution var {"valueS" "frequencyS" "log_size" "log_freq_size"};
append;
close log_size_distribution;

/* for duration */
create log_time_distribution var {"valueT" "frequencyT"};
append;
close log_time_distribution;

/*Data visualization and power-law 1/f distribution*/

title "Frequency of avalanche sizes";
call histogram(Size) scale='PROPORTION';

title "Frequency of avalanche durations (number of iterations before stability)";
call histogram(Time) scale='PROPORTION';

quit;

proc sgplot data=log_size_distribution;
scatter x=valueS y=frequencyS / markerattrs = (symbol=circlefilled);
xaxis type=log label="Avalanche size" grid;
yaxis type=log label="Frequency of size" grid;
title "Log-log distribution of avalanche sizes (1/f noise)";
run;

proc sgplot data=log_time_distribution;
   scatter x=valueT y=frequencyT / markerattrs = (symbol=circlefilled);
   xaxis type=log label="Time" grid;
   yaxis type=log label="Frequency of time" grid;
   title "Log-log distribution of avalanche durations";
run;

/*Regression*/
proc reg data=log_size_distribution(where=(log_Size <= 5));
   model log_freq_size = log_size;
run;
