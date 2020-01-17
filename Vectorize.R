vlog = Vectorize(log)

N = 1e7
x = rexp(N)
tmv = system.time(vlog(x))
tm = system.time(log(x))

tmv/tm
