mutable struct SpinTask
    f::Function
    name::String
    status::Symbol #= :waiting :running :done :error =#
    retstr::Union{String,Nothing}
    output::Union{String,Nothing}
    error::Union{String,Nothing}

    SpinTask(f, name) = new(f, name, :waiting, nothing, nothing, nothing)
end

function _io_like_(out)
    color = get(out, :color, false)
    size = displaysize(out)
    iob = IOBuffer()
    ioc = IOContext(iob, :displaysize=>size, :color=>color)
    return ioc
end

function runtask(t::SpinTask)
    t.status = :running
    output = _io_like_(stdout)

    try
        result = t.f(output)
        if !isnothing(result)
            t.retstr = string(result)
        end
        t.status = :done
    catch e
        bt = catch_backtrace()
        errbuf = _io_like_(stderr)

        if get(errbuf, :color, false)
            printstyled(errbuf, "SUBTASK ERROR: "; bold=true, color=:red)
        else
            print(errbuf, "SUBTASK ERROR: ")
        end

        showerror(errbuf, e)
        Base.show_backtrace(errbuf, bt)
        t.error = String(take!(errbuf.io))
        t.status = :error
    end
    t.output = String(take!(output.io))
    nothing
end

function print_with_newlines(io, s)
    if !isnothing(s) && !isempty(s)
        # println(io)
        print(io, s)
        println(io)
    end
end

run_plain(tasks) = run_plain(stdout, tasks)
function run_plain(io, tasks::Vector{SpinTask})
    pad_names!(tasks)
    printlock = ReentrantLock()
    Threads.@threads for t in tasks
        t.status = :running
        runtask(t)
        @lock printlock begin
            if t.status == :done
                printstyled(io, " ✓ "; color=:green)
                printstyled(io, t.name; bold=true)
                if !isnothing(t.retstr)
                    println(io, " -> ", t.retstr)
                else
                    println(io)
                end
                print_with_newlines(io, t.output)
            elseif t.status == :error
                printstyled(io, " ✗ "; color=:red)
                printstyled(io, t.name; bold=true)
                println(io, " -> ERROR")
                print_with_newlines(io, t.output)
                println(io, t.error)
                println(io)
            end
        end
    end
    if any(t -> t.status == :error, tasks)
        error("One or more tasks failed.")
    end
    nothing
end

const anim_chars = ["◐", "◓", "◑", "◒"]
anim_char(i) = anim_chars[mod1(i, length(anim_chars))]
ansi_moveup(n::Int) = string("\e[", n, "A")
ansi_movecol1 = "\e[1G"
ansi_cleartoend = "\e[0J"
ansi_cleartoendofline = "\e[0K"
ansi_enablecursor = "\e[?25h"
ansi_disablecursor = "\e[?25l"

function pad_names!(tasks::Vector{SpinTask})
    min, max = extrema(length(t.name) for t in tasks)
    if max-min < 5
        # pad names for better appearance
        for t in tasks
            t.name = rpad(t.name, max)
        end
    end
end

function run_fancy(tasks::Vector{SpinTask})
    pad_names!(tasks)

    interrupted = Ref(false)

    printer = Base.errormonitor(Threads.@spawn :interactive begin
        n_tasks = length(tasks)
        t = Timer(0; interval=1/10)
        print(stdout, ansi_disablecursor)
        try
            buf = _io_like_(stdout)
            anim_counter = 0
            while !interrupted[] && any(t -> t.status != :printed, tasks)
                wait(t)
                anim_counter += 1

                print(buf, ansi_cleartoend)

                done = 0
                # print done/error
                for t in tasks
                    if t.status == :printed
                        done += 1
                    elseif t.status == :done
                        printstyled(buf, " ✓ "; color=:green)
                        printstyled(buf, t.name; bold=true)
                        if !isnothing(t.retstr)
                            println(buf, " -> ", t.retstr)
                        else
                            println(buf)
                        end
                        print_with_newlines(buf, t.output)
                        t.status = :printed
                        done += 1
                    elseif t.status == :error
                        printstyled(buf, " ✗ "; color=:red)
                        printstyled(buf, t.name; bold=true)
                        print(buf, " -> ERROR")
                        print_with_newlines(buf, t.output)
                        println(buf, t.error)
                        println(buf)
                        t.status = :printed
                        done +=1
                    end
                end

                # only print progress bar if not all done
                if done < n_tasks
                    # print done/total unicode progress bar
                    width = min(displaysize(stdout)[2], 50)

                    maxlength = length(string(n_tasks))

                    print(buf, " ")
                    progressbar(buf, width - 2*maxlength - 3, done, n_tasks)
                    print(buf, " $(lpad(done, maxlength))/$(n_tasks)")
                    println(buf)

                    # print running tasks
                    running = 0
                    for t in tasks
                        t.status != :running && continue
                        animstep = hash(t.name) + anim_counter
                        print(buf, anim_char(animstep), " ", t.name)
                        println(buf)

                        running += 1
                    end
                end

                s = String(take!(buf.io))
                print(stdout, s)

                if done < n_tasks
                    print(stdout, ansi_moveup(running + 1), ansi_movecol1)
                end

            end
        catch e
            if e isa InterruptException
                interrupted[] = true
            end
            rethrow(e)
        finally
            print(stdout, ansi_enablecursor)
        end
    end)

    Threads.@threads for t in tasks
        if !interrupted[]
            t.status = :running
            try
            runtask(t)
            catch e
                if e isa InterruptException
                    interrupted[] = true
                end
                rethrow(e)
            end
        end
    end

    wait(printer)

    return nothing
end

function progressbar(width, done, total)
    # ioc = _io_like_(stdout)
    ioc = IOContext(IOBuffer())
    progressbar(ioc, width, done, total)
    String(take!(ioc.io))
end
function progressbar(io, width, done, total)
    # ━━━ 0
    # ╸━━ 1
    # ━╺━ 2
    # ━╸━ 3
    # ━━╺ 4
    # ━━━ 5
    steps = width * 2
    step = floor(Int, (done / total) * steps)

    if step == 0
        printstyled(io, "━"^width; color=:light_black)
    elseif step == steps
        printstyled(io, "━"^width; color=:green)
    else
        left_full = div(step,2)
        right_full = width - left_full - 1
        for _ in 1:left_full
            printstyled(io, "━"; color=:green)
        end
        if isodd(step)
            printstyled(io, "╸"; color=:green)
        else
            printstyled(io, "╺"; color=:light_black)
        end
        for _ in 1:right_full
            printstyled(io, "━"; color=:light_black)
        end
    end
    io
end
