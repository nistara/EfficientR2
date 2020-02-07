runTag =
function(args = commandArgs(TRUE))
{
    fns = args[1]
    if(length(args) > 2)
        dir = args[3]
    else
        dir = Sys.getenv("DIR")

    if(dir == "")
        dir = "."
    
    # Compute the suffix.
    suffix = basename(gsub("(functions_|\\.R$)", "", fns))
    print(suffix)
    print(dir)
    source(args[1])
    message("read functions")
    prof = file.path(dir, sprintf("%s.prof", suffix) )
 tm = system.time(source(args[2], verbose = FALSE))    
    print(prof)
    Rprof(prof)
    on.exit(Rprof(NULL))
    tm = system.time(source(args[2], verbose = FALSE))
    Rprof(NULL)
    print(tm)
    p = summaryRprof(prof)
    print(p$by.self)
    save(tm, p, Best, args, file = file.path(dir, sprintf("%s.rda", suffix)))
}
