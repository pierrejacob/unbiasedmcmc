# Here are some examples of calls to varselection.script.R
# arguments are JOB_ID, NRUNS, n, p, SNR, k, m
## WARNING; there are hardcoded paths in varselection.script.R, to be changed by users

# Example 1:
# on 6 cores (-j6), run 10 jobs, with NRUNS=100, n = 500, p = 1000, SNR = 1, k = 0, m = 0
# parallel -j6 'Rscript varselection.script.R {}' ::: {1..10} ::: 100 ::: 500 ::: 1000 ::: 1 ::: 0 ::: 0


# Example 2:
# on 10 cores (-j10), run 10 jobs, with NRUNS=10, n = 500, p = 100, 500 and then 1000, SNR = 3, k = 0, m = 0
# parallel -j10 'Rscript varselection.script.R {}' ::: {1..10} ::: 10 ::: 500 ::: {100,500,1000} ::: 3 ::: 0 ::: 0

