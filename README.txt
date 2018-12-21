How to use IterativeAlignment.py


1) Each matrix you want to work with in a single run must be in the same folder as a ".csv" file. In the given format (see example matrix)


2) Each peptide dataset must be divided into high binders and low binders (or non-binders) and labeled as "high.csv" and "low.csv" respectively.
        a. Similar to the matrices each database that is needed to run in a single run, must be in a folder with each subfolder containing 
           both the high and low files.
3) It is recommended to create a seperate folder for each result to save to, but not necessary.


4) The iterations, cutoff, Dataset Selection, and Gap Penalty all must be positive integers. The perturb percentage must be 0 < x < 1.