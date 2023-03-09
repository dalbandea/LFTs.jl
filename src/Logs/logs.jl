
const global git = Git.git()
const global logdirname = "logs"
const global logfilename = "julia_version_input.txt"
const global gitlogfile = "git_log.txt"

function reproducibility_log(wdir::String = joinpath("./results/trash/", string(Dates.now())))

    logdir = joinpath(wdir, LFTs.logdirname)
    mkpath(logdir)

    logfile = joinpath(logdir, LFTs.logfilename)
    global io_stat = open(logfile, "a")
    InteractiveUtils.versioninfo(io_stat)
    write(io_stat, "\n\n")
    write(io_stat, "$(PROGRAM_FILE) $(join(ARGS, " "))")
    write(io_stat, "\n\n")
    write(io_stat, "RNG_state: ", string(copy(Random.default_rng())))
    write(io_stat, "\n")
    write(io_stat, "Recover with `copy!(Random.default_rng(), state)`")
    close(io_stat)
    cp(PROGRAM_FILE, joinpath(logdir,split(PROGRAM_FILE, "/")[end]))

    # Create patch with non-staged changed (from tracked files)
    write(joinpath(logdir, "gitpatch.patch"), readchomp(`$git diff`))

    global io_stat = open(joinpath(logdir, LFTs.gitlogfile), "w")
    # Print current branch to branchinfo file
    write(io_stat, "Current branch: ")
    write(io_stat, readchomp(`$git rev-parse --abbrev-ref HEAD`))
    write(io_stat, "\n")
    # Print current commit of branch to branchinfo file
    write(io_stat, "Current commit: ")
    write(io_stat, readchomp(`$git rev-parse --short HEAD`))
    write(io_stat, "\n\n")
    # Print current estate of the repository to branchinfo file
    write(io_stat, "Estate of repository:")
    write(io_stat, "\n")
    write(io_stat, readchomp(`$git show-ref`))
    close(io_stat)
end
export reproducibility_log

